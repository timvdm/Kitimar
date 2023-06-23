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

//#define ISOMORPHISM_DFS_RECURSIVE
#define ISOMORPHISM_DFS_ITERATIVE
//#define ISOMORPHISM_DFS_ITERATIVE_OPTIMIZED

#define ISOMORPHISM_MAP_CALLBACK
//#define ISOMORPHISM_MAP_COROUTINE
//#define ISOMORPHISM_MAP_RANGE

#if defined(ISOMORPHISM_DFS_RECURSIVE) + defined(ISOMORPHISM_DFS_ITERATIVE) + defined(ISOMORPHISM_DFS_ITERATIVE_OPTIMIZED) != 1
#error Only one ISOMORPHISM_DFS implementation can be used
#endif

#if defined(ISOMORPHISM_MAP_CALLBACK) + defined(ISOMORPHISM_MAP_COROUTINE) + defined(ISOMORPHISM_MAP_RANGE) != 1
#error Only one ISOMORPHISM_MAP implementation can be used
#endif


// ISOMORPHISM_DFS_BACKTRACK
#if defined(ISOMORPHISM_DFS_RECURSIVE)
    #if defined(ISOMORPHISM_MAP_COROUTINE)
        #define ISOMORPHISM_DFS_BACKTRACK co_return
    #else
        #define ISOMORPHISM_DFS_BACKTRACK return
    #endif
#else
    #define ISOMORPHISM_DFS_BACKTRACK continue
#endif

// ISOMORPHISM_DFS_ABORT
#ifdef ISOMORPHISM_MAP_COROUTINE
    #define ISOMORPHISM_DFS_ABORT co_return
#else
    #define ISOMORPHISM_DFS_ABORT return
#endif

#ifdef ISOMORPHISM_MAP_COROUTINE
#include <cppcoro/recursive_generator.hpp>
#endif

template<auto N>
std::ostream& operator<<(std::ostream &os, const std::array<int, N> &map)
{
    os << "[";
    for (auto i = 0; i < map.size(); ++i)
        os << " " << map[i];
    os << " ]";
    return os;
}

std::ostream& operator<<(std::ostream &os, const std::vector<int> &v)
{
    os << "[ ";
    for (auto i : v)
        os << i << " ";
    os << "]";
    return os;
}

namespace Kitimar::CTSmarts {

    using IsomorphismMapping = std::vector<int>;
    using IsomorphismMappings = std::vector<IsomorphismMapping>;



    enum class MapType {
        Single,
        Unique,
        All
    };

    template<MapType T>
    using MapTypeTag = std::integral_constant<MapType, T>;

    static constexpr auto Single = MapTypeTag<MapType::Single>{};
    static constexpr auto Unique = MapTypeTag<MapType::Unique>{};
    static constexpr auto All    = MapTypeTag<MapType::All>{};

    template<typename SmartsT, MapType Type>
    class Isomorphism
    {

        public:
            using Map = std::array<int, SmartsT::numAtoms>;
            #ifdef ISOMORPHISM_MAP_CALLBACK
            using MapGenerator = IsomorphismMappings;
            using DfsReturnType = void;
            #else
            using MapGenerator = cppcoro::recursive_generator<Map>;
            using DfsReturnType = cppcoro::recursive_generator<Map>;
            #endif

            static constexpr inline auto smarts = SmartsT{};
            static constexpr inline auto dfsBonds = getDfsBonds(smarts);

            static_assert(ctll::size(smarts.bonds));
            static_assert(ctll::size(smarts.bonds) == ctll::size(dfsBonds));

            struct BondIters {
                int begin = 0;
                int bond = -1;
                int end = 0;
            };



            Isomorphism(SmartsT, MapTypeTag<Type>)
            {
                m_degrees = getDegrees<smarts.numAtoms>(smarts.bonds);
                m_map.fill(-1);
            }


            bool match(auto &mol)
            {
                #ifdef ISOMORPHISM_MAP_CALLBACK
                matchDfs(mol, nullptr);
                return isDone();
                #else
                for (const auto &map : matchDfs(mol))
                    return true;
                return false;
                #endif
            }

            auto count(auto &mol, int startAtom = - 1)
            {
                auto n = 0;
                #ifdef ISOMORPHISM_MAP_CALLBACK
                auto cb = [&n] (const auto &array) { ++n; };
                matchDfs(mol, cb, startAtom);
                #else
                for (const auto &map : matchDfs(mol, startAtom))
                    ++n;
                #endif
                return n;
            }

            auto single(auto &mol, int startAtom = -1)
            {
                #ifdef ISOMORPHISM_MAP_CALLBACK
                IsomorphismMapping mapCopy;
                auto cb = [&mapCopy] (const auto &map) {
                    mapCopy.resize(map.size());
                    std::ranges::copy(map, mapCopy.begin());
                };
                matchDfs(mol, cb, startAtom);
                return std::make_tuple(isDone(), mapCopy);
                #else
                IsomorphismMapping mapCopy; // FIXME
                for (const auto &map : matchDfs(mol, startAtom)) {
                    mapCopy.resize(map.size());
                    std::ranges::copy(map, mapCopy.begin());
                    return std::make_tuple(true, mapCopy);
                }
                return std::make_tuple(false, mapCopy);
                #endif
            }

            MapGenerator all(auto &mol, int startAtom = -1)
            {
                #ifdef ISOMORPHISM_MAP_CALLBACK
                IsomorphismMappings maps;
                auto cb = [&maps] (const auto &array) {
                    maps.emplace_back(IsomorphismMapping(array.begin(), array.end()));
                };
                matchDfs(mol, cb, startAtom);
                return maps;
                #else
                co_yield matchDfs(mol,  startAtom);
                #endif
            }

            //
            // Atom
            //

            bool matchAtom(auto &mol, const auto &atom)
            {
                #ifdef ISOMORPHISM_MAP_CALLBACK
                matchDfs(mol, nullptr, get_index(mol, atom));
                return isDone();
                #else
                for (const auto &map : matchDfs(mol, get_index(mol, atom)))
                    return true;
                return false;
                #endif
            }

            bool matchBond(auto &mol, const auto &bond)
            {
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);

                if (matchAtom(mol, source, 0, get<0>(smarts.atoms))) {
                    auto index = get_index(mol, source);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    m_stack[0] = {0, 0, get_degree(mol, source)};
                    #ifdef ISOMORPHISM_MAP_CALLBACK
                    auto callback = [this, &mol, target] (const auto &map) {
                        if (m_map[1] == get_index(mol, target))
                            setDone(true);
                    };
                    matchDfs(mol, callback, index); // FIXME: remove default params
                    if (isDone() && m_map[1] == get_index(mol, target))
                        return true;
                    #else
                    for (const auto &map : matchDfs(mol))
                        if (m_map[1] == get_index(mol, target))
                            return true;
                    #endif
                }

                if (matchAtom(mol, target, 0, get<0>(smarts.atoms))) {
                    auto index = get_index(mol, target);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    m_stack[0] = {0, 0, get_degree(mol, target)};
                    #ifdef ISOMORPHISM_MAP_CALLBACK
                    auto callback = [this, &mol, source] (const auto &map) {
                        if (m_map[1] == get_index(mol, source))
                            setDone(true);
                    };
                    matchDfs(mol, callback, index);
                    return isDone() && m_map[1] == get_index(mol, source);
                    #else
                    for (const auto &map : matchDfs(mol))
                        if (m_map[1] == get_index(mol, source))
                            return true;
                    #endif
                }

                return false;
            }

        private:

            bool matchAtom(auto &mol, const auto &molAtom, int qryAtom, auto atomExpr) const noexcept
            {
                if (get_degree(mol, molAtom) < m_degrees[qryAtom])
                    return false;

                return matchAtomExpr(mol, molAtom, atomExpr);
            }

            constexpr auto makeQueryBondInfo(auto queryBond) const noexcept
            {
                return std::make_tuple(queryBond.source, queryBond.target, queryBond.isCyclic, queryBond.isRingClosure);
            }

            constexpr auto getQueryBondInfo(int queryBondIndex) const noexcept
            {
                switch (queryBondIndex) {
                    case 0:
                        if constexpr (smarts.numBonds > 0)
                            return makeQueryBondInfo(get<0>(dfsBonds));
                    case 1:
                        if constexpr (smarts.numBonds > 1)
                            return makeQueryBondInfo(get<1>(dfsBonds));
                    case 2:
                        if constexpr (smarts.numBonds > 2)
                            return makeQueryBondInfo(get<2>(dfsBonds));
                    case 3:
                        if constexpr (smarts.numBonds > 3)
                            return makeQueryBondInfo(get<3>(dfsBonds));
                    case 4:
                        if constexpr (smarts.numBonds > 4)
                            return makeQueryBondInfo(get<4>(dfsBonds));
                    case 5:
                        if constexpr (smarts.numBonds > 5)
                            return makeQueryBondInfo(get<5>(dfsBonds));
                    case 6:
                        if constexpr (smarts.numBonds > 6)
                            return makeQueryBondInfo(get<6>(dfsBonds));
                    case 7:
                        if constexpr (smarts.numBonds > 7)
                            return makeQueryBondInfo(get<7>(dfsBonds));
                    case 8:
                        if constexpr (smarts.numBonds > 8)
                            return makeQueryBondInfo(get<8>(dfsBonds));
                    case 9:
                        if constexpr (smarts.numBonds > 9)
                            return makeQueryBondInfo(get<9>(dfsBonds));
                    case 10:
                        if constexpr (smarts.numBonds > 10)
                            return makeQueryBondInfo(get<10>(dfsBonds));
                    case 11:
                        if constexpr (smarts.numBonds > 11)
                            return makeQueryBondInfo(get<11>(dfsBonds));
                    case 12:
                        if constexpr (smarts.numBonds > 12)
                            return makeQueryBondInfo(get<12>(dfsBonds));
                    case 13:
                        if constexpr (smarts.numBonds > 13)
                            return makeQueryBondInfo(get<13>(dfsBonds));
                    case 14:
                        if constexpr (smarts.numBonds > 14)
                            return makeQueryBondInfo(get<14>(dfsBonds));
                    case 15:
                        if constexpr (smarts.numBonds > 15)
                            return makeQueryBondInfo(get<15>(dfsBonds));
                    case 16:
                        if constexpr (smarts.numBonds > 16)
                            return makeQueryBondInfo(get<16>(dfsBonds));
                    case 17:
                        if constexpr (smarts.numBonds > 17)
                            return makeQueryBondInfo(get<17>(dfsBonds));
                    case 18:
                        if constexpr (smarts.numBonds > 18)
                            return makeQueryBondInfo(get<18>(dfsBonds));
                    case 19:
                        if constexpr (smarts.numBonds > 19)
                            return makeQueryBondInfo(get<19>(dfsBonds));
                    case 20:
                        if constexpr (smarts.numBonds > 20)
                            return makeQueryBondInfo(get<20>(dfsBonds));
                    default:
                        using R = std::tuple<int, int, bool, bool>;
                        return with_n<ctll::size(dfsBonds), R>(queryBondIndex, [this] (auto i) {
                            auto queryBond = get<i>(dfsBonds);
                            return makeQueryBondInfo(queryBond);
                        });
                }
            }

            constexpr auto getQueryBondSource(int queryBondIndex) const noexcept
            {
                return std::get<0>(getQueryBondInfo(queryBondIndex));
            }

            constexpr auto getQueryBondTarget(int queryBondIndex) const noexcept
            {
                return std::get<1>(getQueryBondInfo(queryBondIndex));
            }

            void debugPoint(int queryBondIndex)
            {
                std::cout << "matchDfs<\"" << smarts.input() << "\">(queryBondIndex = " << queryBondIndex << "):" << std::endl;
                std::cout << "    m_map: " << m_map << std::endl;            
                std::cout << "    m_mapped: [ ";
                for (auto v : m_mapped)
                    std::cout << v << " ";
                std::cout << "]" << std::endl;

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

                std::cout << "    bond iters: ";
                for (const auto &iters : m_stack)
                    if (iters.bond >= 0)
                        std::cout << "[ " << iters.bond << "/" << iters.end << " ]  ";
                    else
                        std::cout << "[ - ]  ";
                std::cout << std::endl;

                if (m_mapped.size() == m_map.size()) {
                    for (auto i = 0; i < m_mapped.size(); ++i)
                        assert(std::ranges::count(m_map, i) == m_mapped[i]);
                    assert(std::ranges::count(m_map, -1) == std::ranges::count(m_mapped, false));
                }
            }



#ifdef ISOMORPHISM_DFS_RECURSIVE

            template<typename Bonds = decltype(dfsBonds)>            
            DfsReturnType matchDfs(auto &mol,
                                   #ifdef ISOMORPHISM_MAP_CALLBACK
                                   auto callback,
                                   #endif
                                   int startAtom = -1, Bonds bonds = dfsBonds)
            {
                if constexpr (!smarts.numBonds)
                    ISOMORPHISM_DFS_ABORT;

                if (isDone())
                    ISOMORPHISM_DFS_ABORT;

                constexpr auto queryBondIndex = smarts.numBonds - ctll::size(bonds);
                if constexpr (ISOMORPHISM_DEBUG)
                    debugPoint(queryBondIndex);


                if constexpr (ctll::empty(bonds)) { // Found mapping?
                    //if (stereoMatches())

                    if constexpr(ISOMORPHISM_DEBUG)
                        std::cout << "    found map: " << m_map << std::endl;

                    if constexpr (Type == MapType::Unique) {
                        // create bit mask of atoms (to ensure uniqueness of mapping)
                        std::vector<bool> atoms(num_atoms(mol));
                        for (auto index : m_map)
                            atoms[index] = true;
                        // add the mapping to the result if it is unique
                        auto hash = std::hash<std::vector<bool>>()(atoms);
                        if (m_maps.find(hash) != m_maps.end())
                            ISOMORPHISM_DFS_BACKTRACK;
                        m_maps.insert(hash);
                    }
                    #ifdef ISOMORPHISM_MAP_CALLBACK
                    addMapping(callback);
                    #else
                    co_yield m_map;
                    #endif
                } else {
                    auto queryBond = ctll::front(bonds);
                    auto querySource = queryBond.source;
                    auto queryTarget = queryBond.target;

                    if constexpr (queryBond.isRingClosure) { // Ring closure?
                        auto source = get_atom(mol, m_map[querySource]);
                        auto target = get_atom(mol, m_map[queryTarget]);
                        if (!m_mapped[get_index(mol, target)])
                            ISOMORPHISM_DFS_BACKTRACK;
                        auto bond = Molecule::get_bond(mol, source, target);
                        if (matchBondExpr(mol, bond, queryBond.bondExpr))
                            #ifdef ISOMORPHISM_MAP_CALLBACK
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #else
                            co_yield matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #endif

                    } else if (m_map[querySource] != -1) { // Source atom mapped?

                        auto source = get_atom(mol, m_map[querySource]);

                        for (auto bond : get_bonds(mol, source)) {

                            if constexpr (ISOMORPHISM_DEBUG)
                                if (m_stack[queryBondIndex].bond)
                                    debugPoint(queryBondIndex);


                            auto target = Molecule::get_nbr(mol, bond, source);
                            auto targetIndex = get_index(mol, target);

                            assert(targetIndex < m_mapped.size());

                            if (m_mapped[targetIndex]) {
                                if constexpr (ISOMORPHISM_DEBUG) {
                                    std::cout << "    target already mapped: " << targetIndex << std::endl;
                                    ++m_stack[queryBondIndex].bond;
                                }
                                continue;
                            }

                            if constexpr (queryBond.isCyclic) {
                                if (!is_cyclic_bond(mol, bond)) {
                                    if constexpr (ISOMORPHISM_DEBUG)
                                        ++m_stack[queryBondIndex].bond;
                                    continue;
                                }
                            }
                            if (!matchBondExpr(mol, bond, queryBond.bondExpr)) {
                                if constexpr (ISOMORPHISM_DEBUG) {
                                    std::cout << "    bond does not match" << std::endl;
                                    ++m_stack[queryBondIndex].bond;
                                }
                                continue;
                            }

                            // match target atom
                            if (!matchAtom(mol, target, queryTarget, queryBond.targetExpr)) {
                                if constexpr (ISOMORPHISM_DEBUG) {
                                    std::cout << "    atom does not match" << std::endl;
                                    ++m_stack[queryBondIndex].bond;
                                }
                                continue;
                            }

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    " << queryTarget << " -> " << targetIndex << '\n';

                            // map target atom
                            m_map[queryTarget] = targetIndex;
                            m_mapped[targetIndex] = true;

                            if constexpr (ISOMORPHISM_DEBUG && queryBondIndex + 1 < smarts.numBonds) {
                                auto nextQueryBond = ctll::front(ctll::pop_front(bonds));
                                m_stack[queryBondIndex + 1] = {0, 0, get_degree(mol, get_atom(mol, nextQueryBond.source))};
                            }

                            #ifdef ISOMORPHISM_MAP_CALLBACK
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #else
                            co_yield matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #endif

                            // exit as soon as possible if only one match is required
                            // (single mapping stored in m_map after returning)
                            if (isDone())
                                ISOMORPHISM_DFS_ABORT;

                            // bracktrack target atom                            
                            if constexpr (!queryBond.isRingClosure) {

                                if constexpr (ISOMORPHISM_DEBUG) {
                                    debugPoint(queryBondIndex);
                                    std::cout << "    backtrack: " << m_map[queryTarget] << std::endl;
                                }
                                m_map[queryTarget] = -1;
                                m_mapped[targetIndex] = false;
                                if constexpr (ISOMORPHISM_DEBUG)
                                    if (queryBondIndex + 1 < m_stack.size())
                                        m_stack[queryBondIndex + 1] = {};
                            }

                            if constexpr (ISOMORPHISM_DEBUG)
                                ++m_stack[queryBondIndex].bond;
                        }



                    } else { // No mapped atoms

                        //if constexpr (ISOMORPHISM_DEBUG)
                        //    std::cout << "no mapped atom" << std::endl;

                        if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                            ISOMORPHISM_DFS_ABORT;

                        reset(mol);

                        for (auto atom : get_atoms(mol)) {
                            if (startAtom != -1)
                                atom = get_atom(mol, startAtom);
                            else if (ISOMORPHISM_DEBUG && get_index(mol, atom))
                                debugPoint(queryBondIndex);



                            auto queryBond = ctll::front(dfsBonds);
                            auto queryAtom = queryBond.source;

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    start atom: " <<  get_index(mol, atom) << std::endl;

                            if (!matchAtom(mol, atom, queryAtom, queryBond.sourceExpr))
                                continue;

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    " << queryAtom << " -> " << get_index(mol, atom) << '\n';

                            // map source atom, recursive dfs, backtrack
                            auto index = get_index(mol, atom);
                            m_map[queryAtom] = index;
                            m_mapped[index] = true;

                            if constexpr (ISOMORPHISM_DEBUG)
                                m_stack[queryBondIndex] = {0, 0, get_degree(mol, atom)};

                            #ifdef ISOMORPHISM_MAP_CALLBACK
                            matchDfs(mol, callback, startAtom, dfsBonds);
                            #else
                            co_yield matchDfs(mol, startAtom, dfsBonds);
                            #endif

                            if (isDone())
                                ISOMORPHISM_DFS_ABORT;

                            if constexpr (ISOMORPHISM_DEBUG) {
                                debugPoint(queryBondIndex);
                                std::cout << "    backtrack: " << m_map[queryAtom] << std::endl;
                            }

                            m_map[queryAtom] = -1;
                            m_mapped[index] = false;

                            if constexpr (ISOMORPHISM_DEBUG)
                                m_stack[queryBondIndex] = {};

                            if (startAtom != -1)
                                ISOMORPHISM_DFS_ABORT;;
                        }

                    }
                }
            }


#endif // ISOMORPHISM_DFS_RECURSIVE

#ifdef ISOMORPHISM_DFS_ITERATIVE

            DfsReturnType matchDfs(auto &mol,
                                   #ifdef ISOMORPHISM_MAP_CALLBACK
                                   auto callback,
                                   #endif
                                   int startAtom = -1)
            {
                if constexpr (!smarts.numBonds)
                    ISOMORPHISM_DFS_ABORT;

                if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                    ISOMORPHISM_DFS_ABORT;

                auto atomIndex = 0;
                auto queryBondIndex = 0;

                while (true) {

                    if (isDone())
                        ISOMORPHISM_DFS_ABORT;

                    if constexpr (ISOMORPHISM_DEBUG)
                        debugPoint(queryBondIndex);

                    if (queryBondIndex == smarts.numBonds) { // Found mapping?
                        //if (stereoMatches())

                        if constexpr(ISOMORPHISM_DEBUG)
                            std::cout << "    found map: " << m_map << std::endl;

                        bool ignore = false;

                        if constexpr (Type == MapType::Unique) {
                            // create bit mask of atoms (to ensure uniqueness of mapping)
                            std::vector<bool> atoms(num_atoms(mol));
                            for (auto index : m_map)
                                atoms[index] = true;
                            // add the mapping to the result if it is unique
                            auto hash = std::hash<std::vector<bool>>()(atoms);
                            if (m_maps.find(hash) != m_maps.end())
                                ignore = true;
                            else
                                m_maps.insert(hash);
                        }
                        if (!ignore) {
                            #if defined(ISOMORPHISM_MAP_CALLBACK)
                            addMapping(callback);
                            #else
                            co_yield m_map;
                            #endif
                        }

                        if (isDone())
                            ISOMORPHISM_DFS_ABORT;

                        --queryBondIndex;
                        ++m_stack[queryBondIndex].bond;
                        continue;

                    }

                    auto [querySource, queryTarget, isCyclic, isRingClosure] = getQueryBondInfo(queryBondIndex);

                    if (isRingClosure) { // Ring closure?

                        auto source = get_atom(mol, m_map[querySource]);
                        auto target = get_atom(mol, m_map[queryTarget]);
                        if (!m_mapped[get_index(mol, target)]) {
                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    target already mapped" << std::endl;
                            ISOMORPHISM_DFS_BACKTRACK;
                        }
                        auto bond = Molecule::get_bond(mol, source, target);

                        if (matchBond(mol, bond, queryBondIndex))
                            ++queryBondIndex;

                        continue;
                    }

                    if (m_map[querySource] != -1) { // Source atom mapped?

                        auto source = get_atom(mol, m_map[querySource]);

                        assert(m_stack.size());
                        auto &bondIters = m_stack[queryBondIndex];
                        if (bondIters.bond == bondIters.end) { // Last incident molecule bond?

                            if (!queryBondIndex && startAtom != -1)
                                ISOMORPHISM_DFS_ABORT;

                            if (m_map[queryTarget] != -1) {
                                if (queryBondIndex + 1 == smarts.numBonds) {
                                    if constexpr (ISOMORPHISM_DEBUG)
                                        std::cout << "    backtrack (final target): " << m_map[queryTarget] << std::endl;
                                    m_mapped[m_map[queryTarget]] = false;
                                    m_map[queryTarget] = -1;
                                }
                            }

                            m_stack[queryBondIndex] = {};

                            if (queryBondIndex) {
                                --queryBondIndex;
                                ++m_stack[queryBondIndex].bond;

                                auto prevQueryTarget = getQueryBondTarget(queryBondIndex);
                                if (m_map[prevQueryTarget] != -1) {
                                    if constexpr (ISOMORPHISM_DEBUG)
                                        std::cout << "    backtrack (target): " << m_map[prevQueryTarget] << std::endl;
                                    m_mapped[m_map[prevQueryTarget]] = false;
                                    m_map[prevQueryTarget] = -1;
                                }
                            } else {
                                if constexpr (ISOMORPHISM_DEBUG)
                                    std::cout << "    backtrack (source): " << m_map[querySource] << std::endl;
                                m_mapped[m_map[querySource]] = false;
                                m_map[querySource] = -1;
                            }

                            continue;
                        }

                        auto bonds = get_bonds(mol, source);
                        auto bondIter = std::begin(bonds);
                        std::advance(bondIter, bondIters.bond);
                        auto bond = *bondIter;

                        auto target = Molecule::get_nbr(mol, bond, source);
                        auto targetIndex = get_index(mol, target);

                        assert(targetIndex < m_mapped.size());
                        if (m_mapped[targetIndex]) {
                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    target already mapped: " << targetIndex << std::endl;
                            ++m_stack[queryBondIndex].bond;
                            continue;
                        }

                        if (isCyclic)
                            if (!is_cyclic_bond(mol, bond)) {
                                ++m_stack[queryBondIndex].bond;
                                continue;
                            }

                        if (!matchBond(mol, bond, queryBondIndex)) {
                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    bond does not match" << std::endl;
                            ++m_stack[queryBondIndex].bond;
                            continue;
                        }

                        // match target atom
                        if (!matchAtom(mol, target, queryTarget)) {
                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    atom does not match" << std::endl;
                            ++m_stack[queryBondIndex].bond;
                            continue;
                        }

                        if constexpr (ISOMORPHISM_DEBUG)
                            std::cout << "    " << queryTarget << " -> " << targetIndex << '\n';

                        // map target atom
                        m_map[queryTarget] = targetIndex;
                        m_mapped[targetIndex] = true;

                        ++queryBondIndex;

                        if (queryBondIndex < smarts.numBonds) {
                            auto nextQuerySource = getQueryBondSource(queryBondIndex);
                            m_stack[queryBondIndex] = {0, 0, get_degree(mol, get_atom(mol, m_map[nextQuerySource]))};
                        }

                        continue;
                    }

                    // No mapped atoms
                    assert(!queryBondIndex);
                    reset(mol);

                    auto queryBond = ctll::front(dfsBonds);
                    auto queryAtom = queryBond.source;

                    if (atomIndex >=  num_atoms(mol))
                        ISOMORPHISM_DFS_ABORT;
                    auto atom = get_atom(mol, startAtom == -1 ? atomIndex++ : startAtom);

                    m_stack[queryBondIndex] = {0, 0, get_degree(mol, atom)};

                    if constexpr (ISOMORPHISM_DEBUG)
                        std::cout << "    start atom: " <<  get_index(mol, atom) << std::endl;

                    if (!matchAtom(mol, atom, queryAtom, queryBond.sourceExpr)) {
                        if (startAtom != -1)
                            ISOMORPHISM_DFS_ABORT;
                        continue;
                    }

                    if constexpr (ISOMORPHISM_DEBUG)
                        std::cout << "    " << queryAtom << " -> " << get_index(mol, atom) << '\n';

                    // map source atom, recursive dfs, backtrack
                    auto index = get_index(mol, atom);
                    m_map[queryAtom] = index;
                    m_mapped[index] = true;

                } // while true
            }

#endif // ISOMORPHISM_DFS_ITERATIVE


            constexpr auto reset(auto &mol) noexcept
            {
                //if constexpr (ISOMORPHISM_DEBUG)
                //    std::cout << "reset()" << std::endl;
                setDone(false);
                m_mapped.clear();
                m_mapped.resize(num_atoms(mol));
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
            constexpr auto addMapping(Callback callback) noexcept
            {
                if constexpr (Type == MapType::Single)
                    setDone(true);
                if constexpr (!std::is_same_v<std::nullptr_t, Callback>)
                    callback(m_map);
            }




            constexpr auto matchAtom(auto &mol, const auto &atom, std::size_t queryAtomIndex) const noexcept
            {
                if (get_degree(mol, atom) < m_degrees[queryAtomIndex])
                    return false;

                switch (queryAtomIndex) {
                    case 0:
                        if constexpr (smarts.numAtoms > 0)
                            return matchAtomExpr(mol, atom, get<0>(smarts.atoms));
                    case 1:
                        if constexpr (smarts.numAtoms > 1)
                            return matchAtomExpr(mol, atom, get<1>(smarts.atoms));
                    case 2:
                        if constexpr (smarts.numAtoms > 2)
                            return matchAtomExpr(mol, atom, get<2>(smarts.atoms));
                    case 3:
                        if constexpr (smarts.numAtoms > 3)
                            return matchAtomExpr(mol, atom, get<3>(smarts.atoms));
                    case 4:
                        if constexpr (smarts.numAtoms > 4)
                            return matchAtomExpr(mol, atom, get<4>(smarts.atoms));
                    case 5:
                        if constexpr (smarts.numAtoms > 5)
                            return matchAtomExpr(mol, atom, get<5>(smarts.atoms));
                    case 6:
                        if constexpr (smarts.numAtoms > 6)
                            return matchAtomExpr(mol, atom, get<6>(smarts.atoms));
                    case 7:
                        if constexpr (smarts.numAtoms > 7)
                            return matchAtomExpr(mol, atom, get<7>(smarts.atoms));
                    case 8:
                        if constexpr (smarts.numAtoms > 8)
                            return matchAtomExpr(mol, atom, get<8>(smarts.atoms));
                    case 9:
                        if constexpr (smarts.numAtoms > 9)
                            return matchAtomExpr(mol, atom, get<9>(smarts.atoms));
                    case 10:
                        if constexpr (smarts.numAtoms > 10)
                            return matchAtomExpr(mol, atom, get<10>(smarts.atoms));
                    case 11:
                        if constexpr (smarts.numAtoms > 11)
                            return matchAtomExpr(mol, atom, get<11>(smarts.atoms));
                    case 12:
                        if constexpr (smarts.numAtoms > 12)
                            return matchAtomExpr(mol, atom, get<12>(smarts.atoms));
                    case 13:
                        if constexpr (smarts.numAtoms > 13)
                            return matchAtomExpr(mol, atom, get<13>(smarts.atoms));
                    case 14:
                        if constexpr (smarts.numAtoms > 14)
                            return matchAtomExpr(mol, atom, get<14>(smarts.atoms));
                    case 15:
                        if constexpr (smarts.numAtoms > 15)
                            return matchAtomExpr(mol, atom, get<15>(smarts.atoms));
                    case 16:
                        if constexpr (smarts.numAtoms > 16)
                            return matchAtomExpr(mol, atom, get<16>(smarts.atoms));
                    case 17:
                        if constexpr (smarts.numAtoms > 17)
                            return matchAtomExpr(mol, atom, get<17>(smarts.atoms));
                    case 18:
                        if constexpr (smarts.numAtoms > 18)
                            return matchAtomExpr(mol, atom, get<18>(smarts.atoms));
                    case 19:
                        if constexpr (smarts.numAtoms > 19)
                            return matchAtomExpr(mol, atom, get<19>(smarts.atoms));
                    case 20:
                        if constexpr (smarts.numAtoms > 20)
                            return matchAtomExpr(mol, atom, get<20>(smarts.atoms));
                    default:
                        return with_n<ctll::size(smarts.atoms), bool>(queryAtomIndex, [&mol, &atom] (auto i) {
                            return matchAtomExpr(mol, atom, get<i>(smarts.atoms));
                        });
                }
            }

            constexpr auto matchBond(auto &mol, const auto &bond, std::size_t queryBondIndex) const noexcept
            {
                switch (queryBondIndex) {
                    case 0:
                        if constexpr (smarts.numBonds > 0)
                            return matchBondExpr(mol, bond, get<0>(dfsBonds).bondExpr);
                    case 1:
                        if constexpr (smarts.numBonds > 1)
                            return matchBondExpr(mol, bond, get<1>(dfsBonds).bondExpr);
                    case 2:
                        if constexpr (smarts.numBonds > 2)
                            return matchBondExpr(mol, bond, get<2>(dfsBonds).bondExpr);
                    case 3:
                        if constexpr (smarts.numBonds > 3)
                            return matchBondExpr(mol, bond, get<3>(dfsBonds).bondExpr);
                    case 4:
                        if constexpr (smarts.numBonds > 4)
                            return matchBondExpr(mol, bond, get<4>(dfsBonds).bondExpr);
                    case 5:
                        if constexpr (smarts.numBonds > 5)
                            return matchBondExpr(mol, bond, get<5>(dfsBonds).bondExpr);
                    case 6:
                        if constexpr (smarts.numBonds > 6)
                            return matchBondExpr(mol, bond, get<6>(dfsBonds).bondExpr);
                    case 7:
                        if constexpr (smarts.numBonds > 7)
                            return matchBondExpr(mol, bond, get<7>(dfsBonds).bondExpr);
                    case 8:
                        if constexpr (smarts.numBonds > 8)
                            return matchBondExpr(mol, bond, get<8>(dfsBonds).bondExpr);
                    case 9:
                        if constexpr (smarts.numBonds > 9)
                            return matchBondExpr(mol, bond, get<9>(dfsBonds).bondExpr);
                    case 10:
                        if constexpr (smarts.numBonds > 10)
                            return matchBondExpr(mol, bond, get<10>(dfsBonds).bondExpr);
                    case 11:
                        if constexpr (smarts.numBonds > 11)
                            return matchBondExpr(mol, bond, get<11>(dfsBonds).bondExpr);
                    case 12:
                        if constexpr (smarts.numBonds > 12)
                            return matchBondExpr(mol, bond, get<12>(dfsBonds).bondExpr);
                    case 13:
                        if constexpr (smarts.numBonds > 13)
                            return matchBondExpr(mol, bond, get<13>(dfsBonds).bondExpr);
                    case 14:
                        if constexpr (smarts.numBonds > 14)
                            return matchBondExpr(mol, bond, get<14>(dfsBonds).bondExpr);
                    case 15:
                        if constexpr (smarts.numBonds > 15)
                            return matchBondExpr(mol, bond, get<15>(dfsBonds).bondExpr);
                    case 16:
                        if constexpr (smarts.numBonds > 16)
                            return matchBondExpr(mol, bond, get<16>(dfsBonds).bondExpr);
                    case 17:
                        if constexpr (smarts.numBonds > 17)
                            return matchBondExpr(mol, bond, get<17>(dfsBonds).bondExpr);
                    case 18:
                        if constexpr (smarts.numBonds > 18)
                            return matchBondExpr(mol, bond, get<18>(dfsBonds).bondExpr);
                    case 19:
                        if constexpr (smarts.numBonds > 19)
                            return matchBondExpr(mol, bond, get<19>(dfsBonds).bondExpr);
                    case 20:
                        if constexpr (smarts.numBonds > 20)
                            return matchBondExpr(mol, bond, get<20>(dfsBonds).bondExpr);
                    default:
                        return with_n<ctll::size(smarts.bonds), bool>(queryBondIndex, [&mol, &bond] (auto i) {
                            return matchBondExpr(mol, bond, get<i>(dfsBonds).bondExpr);
                        });
                }
            }





            std::array<BondIters, smarts.numBonds> m_stack; // FIXME: debug only

            Map m_map; // current mapping: query atom index -> queried atom index
            std::array<uint8_t, smarts.numAtoms> m_degrees; // degree of query atoms
            std::vector<bool> m_mapped; // current mapping: queried atom index -> true if mapped
            std::conditional_t<Type == MapType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };

} // namespace ctsmarts
