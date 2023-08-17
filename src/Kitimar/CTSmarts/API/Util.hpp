#pragma once

#include <Kitimar/Molecule/Molecule.hpp>

#include <vector>
#include <algorithm>
#include <functional>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename Smarts, auto NumCaptures>
        auto captureAtoms(Molecule::Molecule auto &mol, Smarts, bool found, const auto &map,
                          const std::array<int, NumCaptures> &cap)
        {
            using Atom = decltype(get_atom(mol, 0));
            if constexpr (NumCaptures) {
                std::array<Atom, NumCaptures> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < NumCaptures; ++i)
                        atoms[i] = get_atom(mol, map[cap[i]]);
                return atoms;
            } else {
                std::array<Atom, Smarts::numAtoms> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < Smarts::numAtoms; ++i)
                        atoms[i] = get_atom(mol, map[i]); // FIXME: null atoms...
                return atoms;
            }
        }

        auto captureMatchAtoms(Molecule::Molecule auto &mol, auto smarts, bool found, const auto &map, const auto &cap)
        {
            return std::tuple_cat(std::make_tuple(found), captureAtoms(mol, smarts, found, map, cap));
        }

        template<auto N>
        auto copyCapture(Molecule::Molecule auto &mol, const auto &iso, const std::array<int, N> &cap, const auto &caps) noexcept
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto M = N ? N : iso.smarts.numAtoms;
            std::vector<std::array<Atom, M>> v;
            auto r = caps | std::views::transform([&] (const auto &map) {
                return captureAtoms(mol, iso.smarts, true, map, cap);
            });
            std::ranges::copy(r, std::back_inserter(v));
            return v;
        }

        constexpr std::size_t captureHash(auto &mol, const auto &capture)
        {
            auto atoms = std::vector<bool>(num_atoms(mol));
            for (const auto &atom : capture)
                atoms[get_index(mol, atom)] = true;
            return std::hash<std::vector<bool>>()(atoms);
        }

        template<typename Proj = std::identity>
        constexpr auto numInversions(const auto &permutation, Proj proj = {})
        {
            std::size_t n = 0;
            for (auto i = 0; i < permutation.size(); ++i)
                for (auto j = i + 1; j < permutation.size(); ++j)
                    if (std::invoke(proj, permutation[i]) > std::invoke(proj, permutation[j]))
                        ++n;
            return n;
        }

        constexpr bool singleAtomMatch(auto smarts, auto &mol, const auto &atom)
        {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr);
        }

        constexpr bool singleBondMatchHelper(auto smarts, auto &mol, const auto &bond, const auto &source, const auto &target)
        {
            return matchAtomExpr(mol, source, get<0>(smarts.atoms).expr) && matchAtomExpr(mol, target, get<1>(smarts.atoms).expr);
        }

        // 0 -> no match
        // 1 -> source is SMARTS atom 0, target is SMARTS atom 1
        // 2 -> source is SMARTS atom 1, target is SMARTS atom 0
        constexpr int singleBondMatch(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return 1;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return 2;
            return 0;
        }

        constexpr int singleBondCount(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            auto n = 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                ++n;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                ++n;
            return n;
        }

        template<typename Map>
        constexpr auto singleBondMap(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return std::make_tuple(false, Map{});

            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return std::make_tuple(true, Map{get_index(mol, get_source(mol, bond)),
                                                 get_index(mol, get_target(mol, bond))});
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return std::make_tuple(true, Map{get_index(mol, get_target(mol, bond)),
                                                 get_index(mol, get_source(mol, bond))});

            return std::make_tuple(false, Map{});
        }

        template<typename Map>
        constexpr void singleBondMaps(auto smarts, auto &mol, const auto &bond, std::vector<Map> &maps)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                maps.push_back({get_index(mol, get_source(mol, bond)),
                                get_index(mol, get_target(mol, bond))});
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                maps.push_back({get_index(mol, get_target(mol, bond)),
                                get_index(mol, get_source(mol, bond))});
        }

        constexpr auto singleBondCapture(auto smarts, auto &mol, const auto &bond, int singleBondMatchType)
        {
            constexpr auto cap = captureMapping(smarts);
            if constexpr (cap.size() == 1) {
                switch (singleBondMatchType) {
                    case 1:
                        return cap[0] == 0 ? std::make_tuple(true, get_source(mol, bond)) :
                                             std::make_tuple(true, get_target(mol, bond));
                    case 2:
                        return cap[0] == 0 ? std::make_tuple(true, get_target(mol, bond)) :
                                             std::make_tuple(true, get_source(mol, bond));
                    default:
                        return std::make_tuple(false, null_atom(mol));
                }
            } else {
                if (!singleBondMatchType)
                    return std::make_tuple(false, null_atom(mol), null_atom(mol));

                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                if constexpr (cap.size())
                    if (cap[0] > cap[1])
                        std::swap(source, target);

                if (singleBondMatchType == 1)
                    return std::make_tuple(true, source, target);
                return std::make_tuple(true, target, source);
            }
        }

        template<typename AtomMap>
        constexpr void singleBondCaptures(auto smarts, auto &mol, const auto &bond, const auto &cap, std::vector<AtomMap> &maps, bool unique)
        {
            using IndexMap = std::array<decltype(get_index(mol, get_atom(mol, 0))), 2>;
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return;

            if (singleBondMatchHelper(smarts, mol, bond, source, target)) {
                auto map = IndexMap{get_index(mol, get_source(mol, bond)),
                                    get_index(mol, get_target(mol, bond))};
                maps.push_back(impl::captureAtoms(mol, smarts, true, map, cap));
                if (unique)
                    return;
            }

            if (singleBondMatchHelper(smarts, mol, bond, target, source)) {
                auto map = IndexMap{get_index(mol, get_target(mol, bond)),
                                    get_index(mol, get_source(mol, bond))};
                maps.push_back(impl::captureAtoms(mol, smarts, true, map, cap));
            }
        }

        constexpr auto centralAtomMap(auto smarts, auto &mol, const auto &atom)
        {
            std::array<int, smarts.numAtoms> map;
        }

    } // namespace impl

} // mamespace Kitimar::CTSmarts
