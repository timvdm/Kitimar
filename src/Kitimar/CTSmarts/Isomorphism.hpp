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

#define ISOMORPHISM_DFS_RECURSIVE
//#define ISOMORPHISM_DFS_WITH_N
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

// ISOMORPHISM_DFS_RETURN_TYPE
#ifdef ISOMORPHISM_MAP_COROUTINE
    #define ISOMORPHISM_DFS_RETURN_TYPE cppcoro::recursive_generator<Map>
#else
    #define ISOMORPHISM_DFS_RETURN_TYPE void
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
#if defined(ISOMORPHISM_DFS_RECURSIVE)
    #if defined(ISOMORPHISM_MAP_COROUTINE)
        #define ISOMORPHISM_DFS_ABORT co_return
    #else
        #define ISOMORPHISM_DFS_ABORT return
    #endif
#endif

#ifdef ISOMORPHISM_MAP_COROUTINE
#include <cppcoro/recursive_generator.hpp>
#endif


namespace Kitimar::CTSmarts {

    using IsomorphismMapping = std::vector<int>;
    using IsomorphismMappings = std::vector<IsomorphismMapping>;

    template<auto N>
    std::ostream& operator<<(std::ostream &os, const std::array<int, N> &map)
    {
        os << "[";
        for (auto i = 0; i < map.size(); ++i)
            os << " " << map[i];
        os << " ]";
        return os;
    }

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
            #else
            using MapGenerator = cppcoro::recursive_generator<Map>;
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
                    matchDfs(mol, callback, -1, dfsBonds); // FIXME: remove default params
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
                    matchDfs(mol, callback, -1, dfsBonds);
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

            #if defined(ISOMORPHISM_DFS_RECURSIVE)
            template<typename Bonds = decltype(dfsBonds)>
            #endif
            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
            void matchDfs(auto &mol, auto callback, int startAtom = -1, Bonds bonds = dfsBonds)
            #elif defined(ISOMORPHISM_DFS_RECURSIVE) // && !defined(ISOMORPHISM_MAP_CALLBACK)
            ISOMORPHISM_DFS_RETURN_TYPE matchDfs(auto &mol, int startAtom = -1, Bonds bonds = dfsBonds)
            #elif defined(ISOMORPHISM_MAP_CALLBACK) // && !defined(ISOMORPHISM_DFS_RECURSIVE)
            void matchDfs(auto &mol, auto callback, int startAtom = -1)
            #endif
            {
                if constexpr (!smarts.numBonds)
                    ISOMORPHISM_DFS_ABORT;

                if (isDone())
                    ISOMORPHISM_DFS_ABORT;

                if constexpr (ISOMORPHISM_DEBUG)
                    std::cout << "matchDfs()" << std::endl;

                if constexpr (ctll::empty(bonds)) { // Found mapping?
                    //if (stereoMatches())

                    if constexpr(ISOMORPHISM_DEBUG)
                        std::cout << "found map: " << m_map << std::endl;

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
                    #if defined(ISOMORPHISM_MAP_CALLBACK)
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
                            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #elif defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_COROUTINE)
                            co_yield matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #endif

                    } else if (m_map[querySource] != -1) { // Source atom mapped?

                        auto source = get_atom(mol, m_map[querySource]);

                        for (auto bond : get_bonds(mol, source)) {
                            if constexpr (queryBond.isCyclic)
                                if (!is_cyclic_bond(mol, bond))
                                    continue;
                            if (!matchBondExpr(mol, bond, queryBond.bondExpr))
                                continue;

                            auto target = Molecule::get_nbr(mol, bond, source);
                            auto targetIndex = get_index(mol, target);

                            assert(targetIndex < m_mapped.size());
                            if (m_mapped[targetIndex])
                                continue;

                            // match target atom
                            if (!matchAtom(mol, target, queryTarget, queryBond.targetExpr))
                                continue;

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "    " << queryTarget << " -> " << targetIndex << '\n';

                            // map target atom
                            m_map[queryTarget] = targetIndex;
                            m_mapped[targetIndex] = true;

                            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            #elif defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_COROUTINE)
                            co_yield matchDfs(mol, startAtom, ctll::pop_front(bonds));
                            #endif

                            // exit as soon as possible if only one match is required
                            // (single mapping stored in m_map after returning)
                            if (isDone())
                                ISOMORPHISM_DFS_ABORT;

                            // bracktrack target atom
                            if constexpr (!queryBond.isRingClosure) {
                                m_map[queryTarget] = -1;
                                m_mapped[targetIndex] = false;
                            }
                        }

                    } else { // No mapped atoms

                        if constexpr (ISOMORPHISM_DEBUG)
                            std::cout << "no mapped atom" << std::endl;

                        if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                            ISOMORPHISM_DFS_ABORT;

                        reset(mol);

                        for (auto atom : get_atoms(mol)) {
                            if (startAtom != -1)
                                atom = get_atom(mol, startAtom);

                            auto queryBond = ctll::front(dfsBonds);
                            auto queryAtom = queryBond.source;

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << "start atom: " <<  get_index(mol, atom) << std::endl;

                            if (!matchAtom(mol, atom, queryAtom, queryBond.sourceExpr))
                                continue;

                            if constexpr (ISOMORPHISM_DEBUG)
                                std::cout << queryAtom << " -> " << get_index(mol, atom) << '\n';

                            // map source atom, recursive dfs, backtrack
                            auto index = get_index(mol, atom);
                            m_map[queryAtom] = index;
                            m_mapped[index] = true;

                            #if defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_CALLBACK)
                            matchDfs(mol, callback, startAtom, dfsBonds);
                            #elif defined(ISOMORPHISM_DFS_RECURSIVE) && defined(ISOMORPHISM_MAP_COROUTINE)
                            co_yield matchDfs(mol, startAtom, dfsBonds);
                            #endif

                            if (isDone())
                                ISOMORPHISM_DFS_ABORT;
                            m_map[queryAtom] = -1;
                            m_mapped[index] = false;

                            if (startAtom != -1)
                                ISOMORPHISM_DFS_ABORT;;
                        }

                    }
                }
            }

            constexpr auto reset(auto &mol) noexcept
            {
                if constexpr (ISOMORPHISM_DEBUG)
                    std::cout << "reset()" << std::endl;
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



            /*
            constexpr auto matchAtom(auto &mol, const auto &atom, std::size_t queryAtomIndex) const noexcept
            {
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
            */






            Map m_map; // current mapping: query atom index -> queried atom index
            std::array<uint8_t, smarts.numAtoms> m_degrees; // degree of query atoms
            std::vector<bool> m_mapped; // current mapping: queried atom index -> true if mapped
            std::conditional_t<Type == MapType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };

} // namespace ctsmarts
