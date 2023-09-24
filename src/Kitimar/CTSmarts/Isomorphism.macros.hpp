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


#define ISOMORPHISM_DEBUG 1

//#define ISOMORPHISM_DFS_RECURSIVE
#define ISOMORPHISM_DFS_WITH_N
//#define ISOMORPHISM_DFS_OPTIMIZED_WITH_N

#define ISOMORPHISM_MAP_CALLBACK
//#define ISOMORPHISM_MAP_COROUTINE
//#define ISOMORPHISM_MAP_RANGE

#if defined(ISOMORPHISM_DFS_RECURSIVE) + defined(ISOMORPHISM_DFS_WITH_N) + defined(ISOMORPHISM_DFS_OPTIMIZED_WITH_N) != 1
#error Only one ISOMORPHISM_DFS implementation can be used
#endif

#if defined(ISOMORPHISM_MAP_CALLBACK) + defined(ISOMORPHISM_MAP_COROUTINE) + defined(ISOMORPHISM_MAP_RANGE) != 1
#error Only one ISOMORPHISM_MAP implementation can be used
#endif

//// ISOMORPHISM_DFS_NEXT_BOND
//#ifdef ISOMORPHISM_DFS_RECURSIVE
//    #define ISOMORPHISM_DFS_NEXT_BOND continue
//#else
//    #define ISOMORPHISM_DFS_NEXT_BOND continue
//#endif

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
#if defined(ISOMORPHISM_DFS_RECURSIVE)
    #if defined(ISOMORPHISM_MAP_COROUTINE)
        #define ISOMORPHISM_DFS_ABORT co_return
    #else
        #define ISOMORPHISM_DFS_ABORT return
    #endif
#else
    #define ISOMORPHISM_DFS_ABORT return  // FIXME: simplify logic? ifdef ISOMORPHISM_DFS_RECURSIVE else ...
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
                co_yield matchDfs(mol,  get_index(mol, startAtom));
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
                    #ifdef ISOMORPHISM_MAP_CALLBACK
                    auto callback = [this, &mol, target] (const auto &map) {
                        if (m_map[1] == get_index(mol, target))
                            setDone(true);
                    };
                    matchDfs(mol, callback);
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
                    #ifdef ISOMORPHISM_MAP_CALLBACK
                    auto callback = [this, &mol, source] (const auto &map) {
                        if (m_map[1] == get_index(mol, source))
                            setDone(true);
                    };
                    matchDfs(mol, callback);
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

            bool matchAtom(auto &mol, const auto &molAtom, int qryAtom, auto atomExpr) noexcept
            {
                if (get_degree(mol, molAtom) < m_degrees[qryAtom])
                    return false;

                return matchAtomExpr(mol, molAtom, atomExpr);
            }

            void next(auto &mol, auto begin, auto end, auto bonds)
            {
                // - begin not null
                //     - begin == end -> isDone, backtrack, return (null, null, pop(bonds))
                //     - match bond
                //         - match: return (begin+1, end, pop(bonds))
                // - get_source()
                //     - unmapped -> map atom


            }

            constexpr auto queryBondInfo(auto queryBondIndex)
            {
                using R = std::tuple<int, int, bool, bool>;
                return with_n<ctll::size(dfsBonds), R>(queryBondIndex, [] (auto i) {
                    auto queryBond = get<i>(dfsBonds);
                    return std::make_tuple(queryBond.source, queryBond.target, queryBond.isCyclic, queryBond.isRingClosure);
                });
            }

            #ifdef ISOMORPHISM_DFS_RECURSIVE
            template<typename Bonds = decltype(dfsBonds)>
            #endif
            DfsReturnType matchDfs(auto &mol,
            #ifdef ISOMORPHISM_MAP_CALLBACK
                                   auto callback,
            #endif
                                   int startAtom = -1
            #ifdef ISOMORPHISM_DFS_RECURSIVE
                                    , Bonds bonds = dfsBonds
            #endif
            )
            {
                #ifdef ISOMORPHISM_DFS_RECURSIVE
                auto queryBondIndex = ctll::size(dfsBonds) - ctll::size(bonds);
                #else
                if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                    ISOMORPHISM_DFS_ABORT;


                auto ITER = 0;
                auto atomIndex = 0;
                auto queryBondIndex = 0;
                auto depth = 0;





                struct BondIters {
                    int begin = 0;
                    int bond = 0;
                    int end = 0;
                    /*
                    std::ranges::iterator_t<decltype(get_bonds(mol, 0))> begin = {};
                    std::ranges::iterator_t<decltype(get_bonds(mol, 0))> bond = {};
                    std::ranges::sentinel_t<decltype(get_bonds(mol, 0))> end = {};
                    */
                };
                //std::array<BondIters, smarts.numBonds> bondItersStack2 = {};
                std::vector<BondIters> bondItersStack;

                while (true) {
                    //if (atomIndex >= num_atoms(mol))
                    //    ISOMORPHISM_DFS_ABORT;
                    ++ITER;
                    //assert(ITER < 20);
                #endif

                if constexpr (!smarts.numBonds)
                    ISOMORPHISM_DFS_ABORT;

                if (isDone())
                    ISOMORPHISM_DFS_ABORT;

                if constexpr (ISOMORPHISM_DEBUG) {
                    std::cout << "matchDfs(queryBondIndex = " << queryBondIndex << ")" << std::endl;
                    std::cout << "    m_map: " << m_map << std::endl;
                    std::cout << "    bond iters: ";
                    for (const auto &iters : bondItersStack)
                        std::cout << "[ " << iters.bond << "/" << iters.end << " ]  ";
                    std::cout << std::endl;
                    //std::cout << "    detph: " << depth << std::endl;
                }

                // Found mapping?
                #ifdef ISOMORPHISM_DFS_RECURSIVE
                if constexpr (ctll::empty(bonds)) {
                #else
                if (queryBondIndex == smarts.numBonds) {
                #endif
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

                    auto [querySource, queryTarget, isCyclic, isRingClosure] = queryBondInfo(queryBondIndex - 1);


                    std::cout << "    backtrack target: " << m_map[queryTarget] << std::endl;


                    m_mapped[m_map[queryTarget]] = false;
                    m_map[queryTarget] = -1;
                    --queryBondIndex;




                    //if (startAtom != -1)
                    //    return;

                } else {
                    #ifdef ISOMORPHISM_DFS_RECURSIVE
                    auto queryBond = ctll::front(bonds); // constexpr?
                    auto querySource = queryBond.source;
                    auto queryTarget = queryBond.target;
                    auto isCyclic = queryBond.isCyclic;
                    auto isRingClosure = queryBond.isRingClosure;
                    #else
                    auto [querySource, queryTarget, isCyclic, isRingClosure] = queryBondInfo(queryBondIndex);
                    #endif

                    if constexpr(ISOMORPHISM_DEBUG)
                        std::cout << "    queryBond: source = " << querySource << ", target = " << queryTarget
                                  << ", cyclic = " << isCyclic << ", isRingclosure = " << isRingClosure << std::endl;

                    // Ring closure?
                    #ifdef ISOMORPHISM_DFS_RECURSIVE
                    if constexpr (isRingClosure) {
                    #else
                    if (isRingClosure) {
                    #endif
                        auto source = get_atom(mol, m_map[querySource]);
                        auto target = get_atom(mol, m_map[queryTarget]);
                        if (!m_mapped[get_index(mol, target)])
                            ISOMORPHISM_DFS_BACKTRACK;
                        auto bond = Molecule::get_bond(mol, source, target);
                        #ifdef ISOMORPHISM_DFS_RECURSIVE
                        if (matchBondExpr(mol, bond, queryBond.bondExpr)) {
                        #else
                        if (matchBond(mol, bond, queryBondIndex)) {
                        #endif
                            #ifdef ISOMORPHISM_MAP_COROUTINE
                            co_yield
                            #endif
                            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #elif defined(ISOMORPHISM_DFS_RECURSIVE)
                            matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #else
                            ++queryBondIndex;
                            continue;
                            #endif
                        }

                    } else if (m_map[querySource] != -1) { // Source atom mapped?


                        auto source = get_atom(mol, m_map[querySource]);

                        //auto &bondIters = bondItersStack[querySource];
                        assert(bondItersStack.size());
                        auto &bondIters = bondItersStack.back();
                        if constexpr (ISOMORPHISM_DEBUG)
                            std::cout << "    bonds: " << bondIters.bond << " / " << bondIters.end << std::endl;


                        if constexpr (ISOMORPHISM_DEBUG)
                            std::cout << "    mol source: " << get_index(mol, source)  << std::endl;

                        #ifdef ISOMORPHISM_DFS_RECURSIVE
                        for (auto bond : get_bonds(mol, source)) {
                            if constexpr (isCyclic)
                                if (!is_ring_bond(mol, bond))
                                    continue;                            
                            if (!matchBondExpr(mol, bond, queryBond.bondExpr))
                                continue;
                        #else
                            //assert(depth < bondItersStack.size());
                            //auto &bondIters = bondItersStack[depth];
                            if (bondIters.bond == bondIters.end) {
                                //std::cout << "    backtrack target: " << m_map[queryTarget] << std::endl;
                                //m_mapped[m_map[queryTarget]] = false;
                                //m_map[queryTarget] = -1;
                                bondItersStack.pop_back();
                                if (!queryBondIndex && startAtom != -1)
                                    return;
                                if (!bondItersStack.empty())
                                    --queryBondIndex;
                                else {
                                    std::cout << "    backtrack source: " << m_map[querySource] << std::endl;
                                    m_mapped[m_map[querySource]] = false;
                                    m_map[querySource] = -1;
                                }

                                continue;
                                /*
                                if (depth) {
                                    --depth;
                                    std::cout << "    backtrack: " << queryTarget << std::endl;
                                    continue;
                                } else
                                    ISOMORPHISM_DFS_ABORT;
                                */
                            }
                            //auto bond = *bondIters.bond;
                            auto bonds = get_bonds(mol, source);
                            auto bondIter = std::begin(bonds);
                            std::advance(bondIter, bondIters.bond++);
                            auto bond = *bondIter;
                            if (isCyclic)
                                if (!is_ring_bond(mol, bond))
                                    continue;
                        #endif



                            auto target = Molecule::get_nbr(mol, bond, source);
                            auto targetIndex = get_index(mol, target);

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    mol target: " << targetIndex << std::endl;

                            assert(targetIndex < m_mapped.size());
                            if (m_mapped[targetIndex]) {
                                if constexpr (ISOMORPHISM_DEBUG)
                                    std::cout << "    target already mapped" << std::endl;
                                continue;
                            }

                            #ifdef ISOMORPHISM_DFS_RECURSIVE
                            if (!matchBondExpr(mol, bond, queryBond.bondExpr))
                                continue;
                            #else
                            if (!matchBond(mol, bond, queryBondIndex)) {
                                if constexpr (ISOMORPHISM_DEBUG)
                                    std::cout << "    bond does not match" << std::endl;
                                continue;
                            }
                            #endif



                            // match target atom
                            #ifdef ISOMORPHISM_DFS_RECURSIVE
                            if (!matchAtom(mol, target, queryTarget, queryBond.targetExpr))
                                continue;
                            #else
                            if (!matchAtom(mol, target, queryTarget)) {
                                if constexpr (ISOMORPHISM_DEBUG)
                                    std::cout << "    atom does not match" << std::endl;
                                continue;
                            }
                            #endif

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    " << queryTarget << " -> " << targetIndex << '\n';

                            // map target atom
                            m_map[queryTarget] = targetIndex;
                            m_mapped[targetIndex] = true;



                            #ifdef ISOMORPHISM_MAP_COROUTINE
                            co_yield
                            #endif
                            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #elif defined(ISOMORPHISM_DFS_RECURSIVE)
                            matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #else
                            ++queryBondIndex;
                            //if (depth + 1 < bondItersStack.size()) {
                                //auto bonds = get_bonds(mol, target);
                                //assert(queryBondIndex < bondItersStack2.size());

                            if (queryBondIndex < smarts.numBonds) {
                                auto [querySource2, queryTarget2, isCyclic2, isRingClosure2] = queryBondInfo(queryBondIndex);
                                bondItersStack.push_back({0, 0, get_degree(mol, get_atom(mol, m_map[querySource2]))});
                            }

                                //bondItersStack.push_back({0, 0, get_degree(mol, target)});
                                //auto &bondIters2 = bondItersStack2[queryBondIndex];
                                //auto &bondIters2 = bondItersStack[queryTarget];
                                //bondIters2.begin = 0;
                                //bondIters2.bond = 0;
                                //bondIters2.end = get_degree(mol, target);
                                /*
                                bondIters2.begin = std::begin(bonds);
                                bondIters2.bond = std::begin(bonds);
                                bondIters2.end = std::end(bonds);
                                */
                            //}
                            continue;
                            #endif

                            // exit as soon as possible if only one match is required
                            // (single mapping stored in m_map after returning)
                            if (isDone())
                                ISOMORPHISM_DFS_ABORT;

                            // bracktrack target atom
                        #ifdef ISOMORPHISM_DFS_RECURSIVE
                            if constexpr (!isRingClosure) {
                                m_map[queryTarget] = -1;
                                m_mapped[targetIndex] = false;
                            }
                        }
                        #endif

                    } else { // No mapped atoms

                        //assert(!queryBondIndex);

                        //if constexpr (ISOMORPHISM_DEBUG)
                        //    std::cout << "    no mapped atom" << std::endl;

                        #ifdef ISOMORPHISM_DFS_RECURSIVE
                        if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                            ISOMORPHISM_DFS_ABORT;
                        #endif

                        reset(mol);

                        auto queryBond = ctll::front(dfsBonds);
                        auto queryAtom = queryBond.source;

                        #ifdef ISOMORPHISM_DFS_RECURSIVE
                        for (auto atom : get_atoms(mol)) {
                            if (startAtom != -1)
                                atom = get_atom(mol, startAtom);
                        #else
                        //assert(atomIndex < num_atoms(mol));
                        if (atomIndex >=  num_atoms(mol))
                            ISOMORPHISM_DFS_ABORT;
                        auto atom = get_atom(mol, startAtom == -1 ? atomIndex++ : startAtom);

                        //auto bonds = get_bonds(mol, atom);
                        //assert(depth < bondItersStack.size());
                        //auto &bondIters = bondItersStack[queryAtom];
                        bondItersStack.push_back({0, 0, get_degree(mol, atom)});
                        /*
                        assert(bondItersStack.size());
                        auto &bondIters = bondItersStack2[queryBondIndex];
                        bondIters.begin = 0;
                        bondIters.bond = 0;
                        bondIters.end = get_degree(mol, atom);
                        */
                        /*
                        bondIters.begin = std::begin(bonds);
                        bondIters.bond = std::begin(bonds);
                        bondIters.end = std::end(bonds);
                        */

                        #endif



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

                            #ifdef ISOMORPHISM_MAP_COROUTINE
                            co_yield
                            #endif
                            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #elif defined(ISOMORPHISM_DFS_RECURSIVE)
                            matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #else
                            //++queryBondIndex;
                            continue;
                            #endif

                            if (isDone())
                                ISOMORPHISM_DFS_ABORT;


                            #ifdef ISOMORPHISM_DFS_RECURSIVE
                            m_map[queryAtom] = -1;
                            m_mapped[index] = false;
                            #endif

                            if (startAtom != -1)
                                ISOMORPHISM_DFS_ABORT;;
                        #ifdef ISOMORPHISM_DFS_RECURSIVE
                        }
                        #endif

                    }
                }

                #ifndef ISOMORPHISM_DFS_RECURSIVE
                } // while true
                #endif
            }

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

                return with_n<ctll::size(smarts.atoms), bool>(queryAtomIndex, [&mol, &atom] (auto i) {
                    return matchAtomExpr(mol, atom, get<i>(smarts.atoms));
                });
            }

            constexpr auto matchBond(auto &mol, const auto &bond, std::size_t queryBondIndex) const noexcept
            {
                return with_n<ctll::size(smarts.bonds), bool>(queryBondIndex, [&mol, &bond] (auto i) {
                    return matchBondExpr(mol, bond, get<i>(smarts.bonds).expr);
                });
            }







            Map m_map; // current mapping: query atom index -> queried atom index
            std::array<uint8_t, smarts.numAtoms> m_degrees; // degree of query atoms
            std::vector<bool> m_mapped; // current mapping: queried atom index -> true if mapped
            std::conditional_t<Type == MapType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };

} // namespace ctsmarts
