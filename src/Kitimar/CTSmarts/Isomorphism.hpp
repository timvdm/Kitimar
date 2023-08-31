 #pragma once

#include "Smarts.hpp"
#include "Config.hpp"
#include "MatchExpr.hpp"

#include "Filter/NumAtomBondFilter.hpp"
#include "Filter/ElementFilter.hpp"
#include "Filter/FilterPolicy.hpp"

#include <Kitimar/Molecule/Molecule.hpp>

#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <variant>
#include <type_traits>
#include <unordered_set>
#include <cassert>

#define ISOMORPHISM_DEBUG 0








namespace Kitimar::CTSmarts {

    enum class SearchType
    {
        Single,
        Unique,
        All
    };

    template<SearchType T>
    using SearchTypeTag = std::integral_constant<SearchType, T>;

    static constexpr auto Single = SearchTypeTag<SearchType::Single>{};
    static constexpr auto Unique = SearchTypeTag<SearchType::Unique>{};
    static constexpr auto All    = SearchTypeTag<SearchType::All>{};


    template<ctll::fixed_string SMARTS>
    bool requiresExplicitHydrogens() noexcept
    {
        return impl::containsExpr(AliphaticAtom<1>{}, Smarts<SMARTS>::atoms);
    }


    struct NoFilterPolicy : FilterPolicy<> {};

    struct UnconditionalFilterPolicy : FilterPolicy<NumAtomBondFilter> {};

    struct ConditionalFilterPolicy : FilterPolicy<> {};


    template<Molecule::Molecule Mol, typename SmartsT, SearchType SearchT, SeedType SeedT, typename Config = DefaultConfig>
    class Isomorphism
    {

        public:
            using Index = decltype(get_index(std::declval<Mol>(), get_atom(std::declval<Mol>(), 0)));
            using Map = IsomorphismMap<Index, SmartsT::numAtoms>;
            using Maps = IsomorphismMaps<Index, SmartsT::numAtoms>;

            static constexpr inline auto invalidIndex = static_cast<Index>(-1);
            static constexpr inline auto smarts = SmartsT{};
            static constexpr inline auto query = Config::Optimizer::create(smarts, SeedTypeTag<SeedT>{});

            static_assert(smarts.numBonds);
            static_assert(ctll::size(smarts.bonds) == ctll::size(query.bonds));

            constexpr const Map& map() const noexcept
            {
                return m_map.map;
            }

            // match

            bool match(Mol &mol)
            {
                matchDfs(mol, nullptr);
                return isDone();
            }

            bool matchAtom(Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                matchDfs(mol, nullptr, get_index(mol, atom));
                return isDone();
            }

            bool matchBond(Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto sourceCallback = [this, targetIndex] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        setDone(true);
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                if (isDone() && m_map.contains(1, targetIndex))
                    return true;

                auto targetCallback = [this, sourceIndex] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        setDone(true);
                };
                matchDfs(mol, targetCallback, targetIndex);
                return isDone() && m_map.contains(1, sourceIndex);
            }

            // count

            auto count(Mol &mol, int startAtom = -1)
            {
                auto n = 0;
                auto cb = [&n] (const auto &array) { ++n; };
                matchDfs(mol, cb, startAtom);
                return n;
            }

            auto countAtom(Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                return count(mol, get_index(mol, atom));
            }

            auto countBond(Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto n = 0;
                auto sourceCallback = [this, targetIndex, &n] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        ++n;
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                auto targetCallback = [this, sourceIndex, &n] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        ++n;
                };
                matchDfs(mol, targetCallback, targetIndex);
                return n;
            }

            // single

            auto single(Mol &mol, int startAtom = -1)
            {
                matchDfs(mol, nullptr, startAtom);
                return std::make_tuple(isDone(), m_map.map());
            }

            auto singleAtom(Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                return single(mol, get_index(mol, atom));
            }

            auto singleBond(Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto sourceCallback = [this, targetIndex] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        setDone(true);
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                if (isDone() && m_map.contains(1, targetIndex))
                    return std::make_tuple(true, m_map.map());

                auto targetCallback = [this, sourceIndex] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        setDone(true);
                };
                matchDfs(mol, targetCallback, targetIndex);
                if (isDone() && m_map.contains(1, sourceIndex))
                    return std::make_tuple(true, m_map.map());

                return std::make_tuple(false, Map{});
            }

            // all

            Maps all(Mol &mol, int startAtom = -1)
            {
                Maps maps;
                auto cb = [&maps] (const auto &map) {
                    maps.push_back(map);
                };
                matchDfs(mol, cb, startAtom);
                return maps;
            }

            auto allAtom(Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                return all(mol, get_index(mol, atom));
            }

            auto allBond(Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                Maps maps;
                auto sourceCallback = [this, targetIndex, &maps] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        maps.push_back(m_map.map());
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                auto targetCallback = [this, sourceIndex, &maps] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        maps.push_back(m_map.map());
                };
                matchDfs(mol, targetCallback, targetIndex);
                return maps;
            }

        private:


            template<typename Arg, typename ...Args>
            void debug(Arg &&arg, Args &&...args)
            {
                #ifdef KITIMAR_WITH_IOSTREAM
                if constexpr (ISOMORPHISM_DEBUG) {
                    std::cout << std::forward<Arg>(arg);
                    ((std::cout << std::forward<Args>(args)), ...);
                    std::cout << '\n';
                }
                #endif // KITIMAR_WITH_IOSTREAM
            }


            bool matchAtom(Mol &mol, const auto &atom, auto queryAtom) const noexcept
            {
                if (get_degree(mol, atom) < query.degrees[queryAtom.index])
                    return false;

                return matchAtomExpr(mol, atom, queryAtom.expr);
            }

            bool matchBond(Mol &mol, const auto &bond, int queryBondIndex, auto queryBond) const noexcept
            {
                if constexpr (queryBond.isCyclic)
                    if (!is_ring_bond(mol, bond))
                        return false;

                return matchBondExpr(mol, bond, queryBond.expr);
            }



            void addAtom(int queryAtomIndex, auto atomIndex) noexcept
            {
                debug("    ", queryAtomIndex, " -> ", atomIndex);
                assert(static_cast<std::size_t>(queryAtomIndex) < m_map.map().size());
                assert(!m_map.containsQueryAtom(queryAtomIndex));
                m_map.add(queryAtomIndex, atomIndex);
            }

            void removeAtom(int queryAtomIndex, auto atomIndex) noexcept
            {
                debug("    backtrack: ", m_map(queryAtomIndex));
                assert(static_cast<std::size_t>(queryAtomIndex) < m_map.map().size());
                assert(m_map.containsQueryAtom(queryAtomIndex));
                m_map.remove(queryAtomIndex, atomIndex);
            }





            constexpr auto makeQueryBondInfo(auto queryBond) const noexcept
            {
                return std::make_tuple(queryBond.source, queryBond.target, queryBond.isCyclic, queryBond.isRingClosure);
            }

            constexpr auto getQueryBondInfo(int queryBondIndex) const noexcept
            {
                using R = std::tuple<int, int, bool, bool>;
                return with_n<ctll::size(query.bonds), R>(queryBondIndex, [this] (auto i) {
                    auto queryBond = get<i>(query.bonds);
                    return makeQueryBondInfo(queryBond);
                });
            }

            #ifdef KITIMAR_WITH_IOSTREAM
            void debugPoint(int queryBondIndex)
            {
                std::cout << "matchDfs<\"" << smarts.input() << "\">(queryBondIndex = " << queryBondIndex << "):" << std::endl;
                std::cout << "    m_map: " << m_map.map << std::endl;
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

                    auto source = m_map.map[querySource];
                    auto target = m_map.map[queryTarget];

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
            #endif // KITIMAR_WITH_IOSTREAM


            template<typename Bonds = decltype(query.bonds)>
            void matchDfs(Mol &mol, auto callback, int startAtom = -1, Bonds bonds = query.bonds)
            {
                if (isDone())
                    return;

                constexpr auto queryBondIndex = smarts.numBonds - ctll::size(bonds);
                if constexpr (ISOMORPHISM_DEBUG)
                    debugPoint(queryBondIndex);

                if constexpr (ctll::empty(bonds)) { // Found mapping?

                    debug("    found map: ", m_map.map());
                    addMapping(mol, callback);

                } else {
                    auto queryBond = ctll::front(bonds);
                    auto querySource = queryBond.source;
                    auto queryTarget = queryBond.target;

                    if constexpr (queryBond.isRingClosure) { // Ring closure?

                        auto source = get_atom(mol, m_map(querySource.index));
                        auto target = get_atom(mol, m_map(queryTarget.index));
                        if (!m_map.containsAtom(get_index(mol, target)))
                            return;
                        auto bond = Molecule::get_bond(mol, source, target);
                        if (bond == null_bond(mol))
                            return;
                        if (matchBondExpr(mol, bond, queryBond.expr))
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));

                    } else if (m_map.containsQueryAtom(querySource.index)) { // Source atom mapped?

                        auto source = get_atom(mol, m_map(querySource.index));

                        for (auto bond : get_bonds(mol, source)) {

                            auto target = Molecule::get_nbr(mol, bond, source);
                            auto targetIndex = get_index(mol, target);

                            if (m_map.containsAtom(targetIndex)) {
                                debug("    target already mapped: ", targetIndex);
                                continue;
                            }

                            // match bond
                            if (!matchBond(mol, bond, queryBondIndex, queryBond)) {
                                debug("    bond does not match");
                                continue;
                            }

                            // match target atom
                            if (!matchAtom(mol, target, queryTarget)) {
                                debug("    atom does not match");
                                continue;
                            }

                            // map target atom
                            addAtom(queryTarget.index, targetIndex);
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            // exit as soon as possible if only one match is required
                            // (single mapping stored in m_map after returning)
                            if (isDone())
                                return;
                            // bracktrack target atom
                            if constexpr (!queryBond.isRingClosure)
                                removeAtom(queryTarget.index, targetIndex);

                        }

                    } else { // No mapped atoms

                        if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                            return;

                        assert(!isDone());

                        assert(static_cast<std::size_t>(std::ranges::count(m_map.map(), invalidIndex)) == m_map.map().size());
                        //if constexpr (std::is_same_v<MappedPolicy<void>, MappedVector<void>>)
                        //    assert(std::ranges::count(m_mapped, true) == 0);

                        m_map.reset(num_atoms(mol));

                        auto queryBond = ctll::front(query.bonds);
                        auto queryAtom = queryBond.source;

                        for (auto atom : get_atoms(mol)) {
                            if (startAtom != -1)
                                atom = get_atom(mol, startAtom);
                            auto index = get_index(mol, atom);

                            debug("    start atom: ",  index);

                            if (!matchAtom(mol, atom, queryAtom))
                                continue;

                            // map source atom, recursive dfs, backtrack
                            addAtom(queryAtom.index, index);
                            matchDfs(mol, callback, startAtom, query.bonds);
                            if (isDone())
                                return;
                            removeAtom(queryAtom.index, index);

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
                m_map.reset(num_atoms(mol));
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
                if constexpr (SearchT == SearchType::Single)
                    setDone(true);
                if constexpr (SearchT == SearchType::Unique) {
                    // create bit mask of atoms (to ensure uniqueness of mapping)
                    std::vector<bool> atoms(num_atoms(mol));
                    for (auto index : m_map.map())
                        atoms[index] = true;
                    // add the mapping to the result if it is unique
                    auto hash = std::hash<std::vector<bool>>()(atoms);
                    if (m_maps.find(hash) != m_maps.end())
                        return;
                    m_maps.insert(hash);
                }
                if constexpr (!std::is_same_v<std::nullptr_t, Callback>)
                    callback(m_map.map());
            }


            Config::template Map<Index, SmartsT::numAtoms> m_map; // current mapping: query atom index -> queried atom index
            std::conditional_t<SearchT == SearchType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
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

} // namespace Kitimar::CTSmarts
