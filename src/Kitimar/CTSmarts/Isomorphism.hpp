 #pragma once

#include "../Molecule/Molecule.hpp"

#include "Smarts.hpp"
#include "SmartsMatch.hpp"

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>
#include <functional>
#include <variant>
#include <type_traits>
#include <unordered_set>
#include <cassert>


#define ISOMORPHISM_DEBUG 0


template<std::integral I, auto N>
std::ostream& operator<<(std::ostream &os, const std::array<I, N> &map)
{
    os << "[";
    for (auto i = 0; i < map.size(); ++i)
        os << " " << map[i];
    os << " ]";
    return os;
}

template<std::integral I>
std::ostream& operator<<(std::ostream &os, const std::vector<I> &v)
{
    os << "[ ";
    for (auto i : v)
        os << i << " ";
    os << "]";
    return os;
}

namespace Kitimar::CTSmarts {

    enum class MapType
    {
        Single,
        Unique,
        All
    };

    template<std::integral Index, ctll::fixed_string SMARTS>
    using IsomorphismMap = std::array<Index, Smarts<SMARTS>::numAtoms>;

    template<std::integral Index, ctll::fixed_string SMARTS>
    using IsomorphismMaps = std::vector<IsomorphismMap<Index, SMARTS>>;


    template<MapType T>
    using MapTypeTag = std::integral_constant<MapType, T>;

    static constexpr auto Single = MapTypeTag<MapType::Single>{};
    static constexpr auto Unique = MapTypeTag<MapType::Unique>{};
    static constexpr auto All    = MapTypeTag<MapType::All>{};

    template<typename Derived>
    class MappedVector
    {
        public:
            constexpr void resetMapped(auto numAtoms)
            {
                m_mapped.clear();
                m_mapped.resize(numAtoms);
            }

            constexpr bool isMapped(auto atomIndex) const noexcept
            {
                assert(atomIndex < m_mapped.size());
                return m_mapped[atomIndex];
            }

            constexpr void addMapped(auto atomIndex) noexcept
            {
                assert(atomIndex < m_mapped.size());
                assert(!m_mapped[atomIndex]);
                m_mapped[atomIndex] = true;
            }

            constexpr void removeMapped(auto atomIndex) noexcept
            {
                assert(atomIndex < m_mapped.size());
                assert(m_mapped[atomIndex]);
                m_mapped[atomIndex] = false;
            }


        private:
            std::vector<uint8_t> m_mapped;
    };

    template<typename Derived>
    class MappedLookup
    {
        public:
            bool isMapped(auto atomIndex) const noexcept
            {
                const auto &map = static_cast<const Derived*>(this)->map();
                return std::ranges::find(map, atomIndex) != map.end();
            }

            constexpr void resetMapped(auto numAtoms) const noexcept {}
            constexpr void addMapped(auto atomIndex) const noexcept {}
            constexpr void removeMapped(auto atomIndex) const noexcept {}
    };




    template<Molecule::Molecule Mol, typename SmartsT, MapType Type, template<typename> class MappedPolicy = MappedVector>
    class Isomorphism : public MappedPolicy<Isomorphism<Mol, SmartsT, Type>>
    {

        public:
            using Index = decltype(get_index(std::declval<Mol>(), get_atom(std::declval<Mol>(), 0)));
            using Map = IsomorphismMap<Index, SmartsT::smarts>;
            using Maps = IsomorphismMaps<Index, SmartsT::smarts>;

            static constexpr inline auto smarts = SmartsT{};
            static constexpr inline auto dfsBonds = getDfsBonds(smarts);


            static_assert(smarts.numBonds);
            static_assert(ctll::size(smarts.bonds) == ctll::size(dfsBonds));

            Isomorphism()
            {
                m_degrees = getDegrees<smarts.numAtoms>(smarts.bonds);
                m_map.fill(-1);
            }

            constexpr const Map& map() const noexcept
            {
                return m_map;
            }

            bool match(Mol &mol)
            {                
                matchDfs(mol, nullptr);
                return isDone();
            }

            auto count(Mol &mol, int startAtom = - 1)
            {
                auto n = 0;                
                auto cb = [&n] (const auto &array) { ++n; };
                matchDfs(mol, cb, startAtom);
                return n;
            }


            auto single(Mol &mol, int startAtom = -1)
            {                
                matchDfs(mol, nullptr, startAtom);
                return std::make_tuple(isDone(), m_map);                
            }

            Maps all(Mol &mol, int startAtom = -1)
            {                
                Maps maps;
                auto cb = [&maps] (const auto &map) {
                    maps.push_back(map);
                };
                matchDfs(mol, cb, startAtom);
                return maps;             
            }

            //
            // Atom
            //

            bool matchAtom(Mol &mol, const auto &atom)
            {                
                matchDfs(mol, nullptr, get_index(mol, atom));
                return isDone();                
            }

            bool matchBond(Mol &mol, const auto &bond)
            {
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto sourceCallback = [this, targetIndex] (const auto &map) {
                    if (m_map[1] == targetIndex)
                        setDone(true);
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                if (isDone() && m_map[1] == targetIndex)
                    return true;

                auto targetCallback = [this, sourceIndex] (const auto &map) {
                    if (m_map[1] == sourceIndex)
                        setDone(true);
                };
                matchDfs(mol, targetCallback, targetIndex);
                return isDone() && m_map[1] == sourceIndex;
            }

        private:

            template<typename Arg, typename ...Args>
            void debug(Arg &&arg, Args &&...args)
            {
                if constexpr (ISOMORPHISM_DEBUG) {
                    std::cout << std::forward<Arg>(arg);
                    ((std::cout << std::forward<Args>(args)), ...);
                    std::cout << '\n';
                }
            }


            bool matchAtom(Mol &mol, const auto &atom, int queryAtomIndex, auto atomExpr) const noexcept
            {
                if (get_degree(mol, atom) < m_degrees[queryAtomIndex])
                    return false;

                return matchAtomExpr(mol, atom, atomExpr);
            }

            bool matchBond(Mol &mol, const auto &bond, int queryBondIndex, auto queryBond) const noexcept
            {
                if constexpr (queryBond.isCyclic)
                    if (!is_cyclic_bond(mol, bond))
                        return false;

                return matchBondExpr(mol, bond, queryBond.bondExpr);
            }



            void addAtom(auto atomIndex, int queryAtomIndex) noexcept
            {
                debug("    ", queryAtomIndex, " -> ", atomIndex);

                assert(queryAtomIndex < m_map.size());
                assert(m_map[queryAtomIndex] == -1);
                m_map[queryAtomIndex] = atomIndex;

                this->addMapped(atomIndex);
            }

            void removeAtom(auto atomIndex, int queryAtomIndex) noexcept
            {
                debug("    backtrack: ", m_map[queryAtomIndex]);

                assert(queryAtomIndex < m_map.size());
                assert(m_map[queryAtomIndex] == atomIndex);
                m_map[queryAtomIndex] = -1;

                this->removeMapped(atomIndex);
            }





            constexpr auto makeQueryBondInfo(auto queryBond) const noexcept
            {
                return std::make_tuple(queryBond.source, queryBond.target, queryBond.isCyclic, queryBond.isRingClosure);
            }

            constexpr auto getQueryBondInfo(int queryBondIndex) const noexcept
            {
                using R = std::tuple<int, int, bool, bool>;
                return with_n<ctll::size(dfsBonds), R>(queryBondIndex, [this] (auto i) {
                    auto queryBond = get<i>(dfsBonds);
                    return makeQueryBondInfo(queryBond);
                });
            }


            void debugPoint(int queryBondIndex)
            {
                std::cout << "matchDfs<\"" << smarts.input() << "\">(queryBondIndex = " << queryBondIndex << "):" << std::endl;
                std::cout << "    m_map: " << m_map << std::endl;
                /*
                std::cout << "    m_mapped: [ ";
                for (auto v : m_mapped)
                    std::cout << v << " ";
                std::cout << "]" << std::endl;
                */

                if (queryBondIndex >= 0 && queryBondIndex < smarts.numBonds) {
                    // query
                    auto [querySource, queryTarget, isCyclic, isRingClosure] = getQueryBondInfo(queryBondIndex);
                    std::cout << "    query bond: " << querySource << " - " << queryTarget << std::endl;

                    auto source = m_map[querySource];
                    auto target = m_map[queryTarget];

                    std::cout << "    mol bond:   ";
                    if (source < 0)
                        std::cout << "? - ";
                    else
                        std::cout << source << " - ";
                    if (target < 0)
                        std::cout << "?" << std::endl;
                    else
                        std::cout << target << std::endl;
                }
                /*
                if (m_mapped.size() == m_map.size()) {
                    for (auto i = 0; i < m_mapped.size(); ++i)
                        assert(std::ranges::count(m_map, i) == m_mapped[i]);
                    assert(std::ranges::count(m_map, -1) == std::ranges::count(m_mapped, false));
                }
                */
            }



            template<typename Bonds = decltype(dfsBonds)>
            void matchDfs(Mol &mol, auto callback, int startAtom = -1, Bonds bonds = dfsBonds)
            {
                if (isDone())
                    return;

                constexpr auto queryBondIndex = smarts.numBonds - ctll::size(bonds);
                if constexpr (ISOMORPHISM_DEBUG)
                    debugPoint(queryBondIndex);

                if constexpr (ctll::empty(bonds)) { // Found mapping?

                    debug("    found map: ", m_map);
                    addMapping(mol, callback);

                } else {
                    auto queryBond = ctll::front(bonds);
                    auto querySource = queryBond.source;
                    auto queryTarget = queryBond.target;

                    if constexpr (queryBond.isRingClosure) { // Ring closure?

                        auto source = get_atom(mol, m_map[querySource]);
                        auto target = get_atom(mol, m_map[queryTarget]);
                        if (!this->isMapped(get_index(mol, target)))
                            return;
                        auto bond = Molecule::get_bond(mol, source, target);
                        if (bond == null_bond(mol))
                            return;
                        if (matchBondExpr(mol, bond, queryBond.bondExpr))                            
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));                            

                    } else if (m_map[querySource] != -1) { // Source atom mapped?

                        auto source = get_atom(mol, m_map[querySource]);

                        for (auto bond : get_bonds(mol, source)) {

                            auto target = Molecule::get_nbr(mol, bond, source);
                            auto targetIndex = get_index(mol, target);

                            if (this->isMapped(targetIndex)) {
                                debug("    target already mapped: ", targetIndex);
                                continue;
                            }

                            // match bond
                            if (!matchBond(mol, bond, queryBondIndex, queryBond)) {
                                debug("    bond does not match");
                                continue;
                            }

                            // match target atom
                            if (!matchAtom(mol, target, queryTarget, queryBond.targetExpr)) {                                
                                debug("    atom does not match");
                                continue;
                            }

                            // map target atom                            
                            addAtom(targetIndex, queryTarget);
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            // exit as soon as possible if only one match is required
                            // (single mapping stored in m_map after returning)
                            if (isDone())
                                return;
                            // bracktrack target atom
                            if constexpr (!queryBond.isRingClosure)
                                removeAtom(targetIndex, queryTarget);

                        }

                    } else { // No mapped atoms

                        if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                            return;

                        assert(!isDone());

                        assert(std::ranges::count(m_map, -1) == m_map.size());
                        //if constexpr (std::is_same_v<MappedPolicy<void>, MappedVector<void>>)
                        //    assert(std::ranges::count(m_mapped, true) == 0);

                        this->resetMapped(num_atoms(mol));

                        for (auto atom : get_atoms(mol)) {
                            if (startAtom != -1)
                                atom = get_atom(mol, startAtom);
                            auto index = get_index(mol, atom);

                            auto queryBond = ctll::front(dfsBonds);
                            auto queryAtom = queryBond.source;

                            debug("    start atom: ",  index);
                            
                            if (!matchAtom(mol, atom, queryAtom, queryBond.sourceExpr))
                                continue;

                            // map source atom, recursive dfs, backtrack
                            addAtom(index, queryAtom);
                            matchDfs(mol, callback, startAtom, dfsBonds);
                            if (isDone())
                                return;
                            removeAtom(index, queryAtom);

                            if (startAtom != -1)
                                return;;
                        }

                    }
                }
            }

            constexpr auto reset(Mol &mol) noexcept
            {                
                //debug("reset()");
                setDone(false);
                this->resetMapped(num_atoms(mol));
            }

            constexpr auto isDone() const noexcept
            {
                return m_done;
            }

            constexpr auto setDone(bool done) noexcept
            {
                m_done = done;
            }

            template<typename Callback>
            constexpr auto addMapping(Molecule::Molecule auto &mol, Callback callback) noexcept
            {
                if constexpr (Type == MapType::Single)
                    setDone(true);
                if constexpr (Type == MapType::Unique) {
                    // create bit mask of atoms (to ensure uniqueness of mapping)
                    std::vector<bool> atoms(num_atoms(mol));
                    for (auto index : m_map)
                        atoms[index] = true;
                    // add the mapping to the result if it is unique
                    auto hash = std::hash<std::vector<bool>>()(atoms);
                    if (m_maps.find(hash) != m_maps.end())
                        return;
                    m_maps.insert(hash);
                }
                if constexpr (!std::is_same_v<std::nullptr_t, Callback>)
                    callback(m_map);
            }


            Map m_map; // current mapping: query atom index -> queried atom index
            std::array<uint8_t, smarts.numAtoms> m_degrees; // degree of query atoms
            std::conditional_t<Type == MapType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };


    /*
    template<Molecule::Molecule Mol, typename SmartsT, MapType Type>
    class CentralAtomIsomorphism
    {
        public:
            using Map = IsomorphismMap<SmartsT::smarts>;
            using Maps = IsomorphismMaps<SmartsT::smarts>;

            static constexpr inline auto smarts = SmartsT{};

            static_assert(smarts.numBonds);

            CentralAtomIsomorphism()
            {
                m_map.fill(-1);
            }

            auto count(Mol &mol, int startAtom = - 1)
            {
                auto n = 0;
                auto cb = [&n] (const auto &array) { ++n; };
                matchCentalAtom(mol, cb, startAtom);
                return n;
            }



        private:
            template<typename Bonds = decltype(smarts.bonds)>
            void matchCentralAtom(Mol &mol, auto callback, int startAtom = -1, Bonds bonds = smarts.bonds)
            {
                auto queryBond = ctll::front(bonds);
                auto querySource = queryBond.source;
                auto queryTarget = queryBond.target;

                //if constexpr (queryTarget == smarts.centralAtom)

                for (auto atom : get_atoms(mol)) {
                    if (startAtom != -1)
                        atom = get_atom(mol, startAtom);
                    auto index = get_index(mol, atom);

                    debug("    start atom: ",  index);

                    if (matchAtom(mol, atom, querySource, queryBond.sourceExpr))
                        continue;

                    if (startAtom != -1)
                        return;;
                }
            }


            Map m_map; // current mapping: query atom index -> queried atom index
    };
    */

} // namespace ctsmarts
