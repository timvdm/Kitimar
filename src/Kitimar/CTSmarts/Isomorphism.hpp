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


#define DEBUG_ISOMORPHISM 0

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
            static constexpr inline auto smarts = SmartsT{};
            static constexpr inline auto dfsBonds = getDfsBonds(smarts);

            static_assert(ctll::size(smarts.bonds) == ctll::size(dfsBonds));

            Isomorphism(SmartsT, MapTypeTag<Type>)
            {
                m_degrees = getDegrees<smarts.numAtoms>(smarts.bonds);
                m_map.fill(-1);
            }

            //
            // Molecule
            //

            bool match(auto &mol)
            {
                //if (num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                if (num_atoms(mol) < smarts.numAtoms)
                    return false;
                matchComponents(mol, nullptr);
                return isDone();
            }

            auto count(auto &mol)
            {                
                auto n = 0;
                matchComponents(mol, [&n] (const auto &array) { ++n; });
                return n;
            }

            auto single(auto &mol)
            {
                IsomorphismMapping map;
                matchComponents(mol, [&map] (const auto &array) {
                    map.resize(array.size());
                    std::ranges::copy(array, map.begin());
                });
                return map;
            }

            auto all(auto &mol)
            {
                IsomorphismMappings maps;
                matchComponents(mol, [&maps] (const auto &array) {
                    maps.emplace_back(IsomorphismMapping(array.begin(), array.end()));

                });
                return maps;
            }

            //
            // Atom
            //

            bool matchAtom(auto &mol, const auto &atom)
            {
                reset(mol);
                matchComponent(mol, atom,  nullptr);
                return isDone();
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
                    auto callback = [this, &mol, target] (const auto &map) {
                        if (m_map[1] == get_index(mol, target))
                            setDone(true);
                    };
                    matchDfs(mol, callback, dfsBonds);
                    if (isDone() && m_map[1] == get_index(mol, target))
                        return true;
                }

                reset(mol); // FIXME: needed?
                if (matchAtom(mol, target, 0, get<0>(smarts.atoms))) {
                    auto index = get_index(mol, target);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    auto callback = [this, &mol, source] (const auto &map) {
                        if (m_map[1] == get_index(mol, source))
                            setDone(true);
                    };
                    matchDfs(mol, callback, dfsBonds);
                    return isDone() && m_map[1] == get_index(mol, source);
                }

                return false;
            }

            auto count(auto &mol, const auto &atom)
            {
                reset(mol);
                auto n = 0;
                matchComponent(mol, atom, [&n] (const auto &array) { ++n; });
                return n;
            }

            auto single(auto &mol, const auto &atom)
            {
                reset(mol);
                IsomorphismMapping map;
                matchComponent(mol, atom, [&map] (const auto &array) {
                    map.resize(array.size());
                    std::ranges::copy(map, map.begin());
                });
                return map;
            }

            auto all(auto &mol, const auto &atom)
            {
                reset(mol);
                IsomorphismMappings maps;
                matchComponent(mol, atom, [&maps] (const auto &array) {
                    maps.emplace_back(IsomorphismMapping(array.begin(), array.end()));

                });
                return maps;
            }


        private:

            bool matchAtom(auto &mol, const auto &molAtom, int qryAtom, auto atomExpr) noexcept
            {
                if (get_degree(mol, molAtom) < m_degrees[qryAtom])
                    return false;

                return matchAtomExpr(mol, molAtom, atomExpr);
            }

            void matchDfs(auto &mol, auto callback, auto bonds)
            {
                static_assert(!ctll::empty(bonds));
                auto queryBond = ctll::front(bonds);
                auto querySource = queryBond.source;
                auto queryTarget = queryBond.target;

                auto source = get_atom(mol, m_map[querySource]);

                //if (DEBUG_ISOMORPHISM)
                //    std::cout << "mapping bond: " << get_index(m_query, querySource) << "-" << get_index(m_query, queryTarget) << std::endl;

                for (auto bond : get_bonds(mol, source)) {

                    if constexpr (queryBond.isCyclic)
                        if (!is_cyclic_bond(mol, bond))
                            continue;

                    if (!matchBondExpr(mol, bond, queryBond.bondExpr))
                        continue;

                    auto target = Molecule::get_nbr(mol, bond, source);
                    auto targetIndex = get_index(mol, target);


                    if constexpr (queryBond.isRingClosure) {
                        if (m_map[queryTarget] != targetIndex)
                            continue;
                    } else {
                        assert(targetIndex < m_mapped.size());
                        if (m_mapped[targetIndex])
                            continue;

                        // match target atom
                        //if (!matchAtomExpr(mol, target, queryBond.targetExpr))
                        if (!matchAtom(mol, target, queryTarget, queryBond.targetExpr))
                            continue;

                        if constexpr (DEBUG_ISOMORPHISM)
                            std::cout << "    " << queryTarget << " -> " << targetIndex << '\n';

                        // map target atom
                        m_map[queryTarget] = targetIndex;
                        m_mapped[targetIndex] = true;
                    }

                    if constexpr (ctll::size(bonds) == 1) {
                        // Found mapping...

                        //if (stereoMatches())

                        //if (DEBUG_ISOMORPHISM)
                        //    std::cout << "found map: " << m_map << std::endl;

                        if constexpr (Type == MapType::Unique) {
                            // create bit mask of atoms (to ensure uniqueness of mapping)
                            std::vector<bool> atoms(num_atoms(mol));
                            for (auto index : m_map)
                                atoms[index] = true;
                            // add the mapping to the result if it is unique
                            auto hash = std::hash<std::vector<bool>>()(atoms);
                            if (m_maps.find(hash) == m_maps.end()) {
                                m_maps.insert(hash);
                                addMapping(callback);
                            }
                        } else
                            addMapping(callback);
                    } else {
                        matchDfs(mol, callback, ctll::pop_front(bonds));
                    }

                    // exit as soon as possible if only one match is required
                    // (single mapping stored in m_map after returning)
                    if (isDone())
                        return;

                    // bracktrack target atom
                    if constexpr (!queryBond.isRingClosure) {
                        m_map[queryTarget] = -1;
                        m_mapped[targetIndex] = false;
                    }
                }
            }


            void matchComponent(auto &mol, const auto &atom, auto callback)
            {
                if constexpr (!smarts.numAtoms)
                    return;

                if (isDone())
                    return;

                if constexpr (!smarts.numBonds) {
                    constexpr auto queryAtom = 0;

                    //if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms)))
                    if (!matchAtom(mol, atom, 0, get<0>(smarts.atoms)))
                        return;

                    if constexpr (DEBUG_ISOMORPHISM)
                        std::cout << queryAtom << " -> " << get_index(mol, atom) << '\n';

                    // map source atom, recursive dfs, backtrack
                    auto index = get_index(mol, atom);
                    m_map[0] = index;
                    m_mapped[index] = true;
                    addMapping(callback);
                    m_map[0] = -1;
                    m_mapped[index] = false;

                } else {
                    auto queryBond = ctll::front(dfsBonds);
                    auto queryAtom = queryBond.source;


                    //if (!matchAtomExpr(mol, atom, queryBond.sourceExpr))
                    if (!matchAtom(mol, atom, queryAtom, queryBond.sourceExpr))
                        return;

                    if constexpr (DEBUG_ISOMORPHISM)
                        std::cout << queryAtom << " -> " << get_index(mol, atom) << '\n';

                    // map source atom, recursive dfs, backtrack
                    auto index = get_index(mol, atom);
                    m_map[queryAtom] = index;
                    m_mapped[index] = true;
                    matchDfs(mol, callback, dfsBonds);
                    m_map[queryAtom] = -1;
                    m_mapped[index] = false;
                }
            }

            void matchComponents(auto &mol, auto callback)
            {
                if constexpr (!smarts.numAtoms)
                    return;

                reset(mol);

                if constexpr (!smarts.numBonds) {
                    for (auto atom : get_atoms(mol))
                        matchComponent(mol, atom, callback);
                } else {
                    // try to match each atom in the molecule against the first atom
                    // epxression in the SMARTS
                    for (auto atom : get_atoms(mol))
                        matchComponent(mol, atom, callback);
                }
            }

            constexpr auto reset(auto &mol) noexcept
            {
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






            std::array<int, smarts.numAtoms> m_map; // current mapping: query atom index -> queried atom index
            std::array<uint8_t, smarts.numAtoms> m_degrees; // degree of query atoms

            //std::array<bool, smarts.numAtoms> m_mapped; // current mapping: queried atom index -> true if mapped
            std::vector<bool> m_mapped; // current mapping: queried atom index -> true if mapped
            std::conditional_t<Type == MapType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };

} // namespace ctsmarts
