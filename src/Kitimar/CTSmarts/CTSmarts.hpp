#pragma once

#include "Isomorphism.hpp"
#include "../Molecule/MockMolecule.hpp" // for testing single header using compiler explorer

namespace Kitimar::CTSmarts {

    namespace detail {

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

        constexpr bool singleAtomMatch(auto smarts, auto &mol, const auto &atom)
        {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms));
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
            if (matchAtomExpr(mol, source, get<0>(smarts.atoms)) && matchAtomExpr(mol, target, get<1>(smarts.atoms)))
                return 1;
            if (matchAtomExpr(mol, source, get<1>(smarts.atoms)) && matchAtomExpr(mol, target, get<0>(smarts.atoms)))
                return 2;
            return 0;
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

        constexpr auto centralAtomMap(auto smarts, auto &mol, const auto &atom)
        {
            std::array<int, smarts.numAtoms> map;
        }


    } // namespace detail


    //
    // CTSmarts::contains<"SMARTS">(mol) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool contains(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol)) // FIXME: use std::ranges::find_if -> check assembly?
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return true;
            return false;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol))
                if (detail::singleBondMatch(smarts, mol, bond))
                    return true;
            return false;
        } else {
            auto iso = Isomorphism{smarts, Single};
            return iso.match(mol);
        }
    }

    //
    // CTSmarts::atom<"SMARTS">(mol, atom) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool atom(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            return detail::singleAtomMatch(smarts, mol, atom);
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms)))
                return false;
            for (auto bond : get_bonds(mol, atom)) {
                if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                    continue;
                if (matchAtomExpr(mol, get_nbr(mol, bond, atom), get<1>(smarts.atoms)))
                    return true;
            }
            return false;
        } else {
            auto iso = Isomorphism{smarts, Single};
            return iso.matchAtom(mol, atom);
        }
    }

    //
    // CTSmarts::bond<"SMARTS">(mol, bond) -> bool
    //

    template<ctll::fixed_string SMARTS>
    constexpr bool bond(Molecule::Molecule auto &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        //std::cout << "CTSmarts::bond<" << smarts.input() << ">(mol, " << get_index(mol, bond) << ")" << std::endl;
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            return detail::singleBondMatch(smarts, mol, bond);
        } else {
            auto iso = Isomorphism{smarts, All};
            return iso.matchBond(mol, bond);
            /*
            for (const auto &map : iso.all(mol)) {
                if (map[queryBond.source] == sourceIndex && map[queryBond.target] == targetIndex)
                    return true;
                if (map[queryBond.source] == targetIndex && map[queryBond.target] == sourceIndex)
                    return true;
            }
            return false;
            */
        }
    }

    //
    // CTSmarts::count<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::integeral
    //

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    constexpr auto count(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        static_assert(M != MapType::Single, "Use CTSmarts::contains<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            auto n = 0;
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    ++n;
            return n;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            auto n = 0;
            for (auto bond : get_bonds(mol))
                if (detail::singleBondMatch(smarts, mol, bond))
                    ++n;
            return n;
        } else {
            auto iso = Isomorphism{smarts, mapType};
            return iso.count(mol);
        }
    }

    //
    // CTSmarts::single<"SMARTS">(mol) -> std::vector<int>
    //

    template<ctll::fixed_string SMARTS>
    constexpr auto single(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return IsomorphismMapping{1, atom};
            return IsomorphismMapping{};
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                switch (detail::singleBondMatch(smarts, mol, bond)) {
                    case 1:
                        return IsomorphismMapping{1, get_source(mol, bond), get_target(mol, bond)};
                    case 2:
                        return IsomorphismMapping{1, get_target(mol, bond), get_source(mol, bond)};
                    default:
                        break;
                }
            }
            return IsomorphismMapping{};
        } else {
            auto iso = Isomorphism{smarts, Single};
            return iso.single(mol);
        }
    }

    //
    // CTSmarts::single<"SMARTS">(mol, atom) -> std::vector<int>
    //

    template<ctll::fixed_string SMARTS>
    constexpr auto single(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        return iso.single(mol, atom);
    }

    //
    // CTSmarts::multi<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::vector<std::vector<int>>
    //

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    constexpr auto multi(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        return iso.all(mol);
    }

    //
    // CTSmarts::capture<"SMARTS">(mol) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS>
    auto capture(Molecule::Molecule auto &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, atom);
            return std::make_tuple(false, null_atom(mol));
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            constexpr auto cap = captureMapping(smarts);
            for (auto bond : get_bonds(mol)) {
                auto matchType = detail::singleBondMatch(smarts, mol, bond);
                if (matchType)
                    return detail::singleBondCapture(smarts, mol, bond, matchType);
            }
            return detail::singleBondCapture(smarts, mol, null_bond(mol), 0);
        } else {
            auto iso = Isomorphism{smarts, Single};
            constexpr auto cap = captureMapping(smarts);
            auto [found, map] = iso.single(mol);
            return detail::captureMatchAtoms(mol, smarts, found, map, cap);
        }
    }

    //
    // CTSmarts::captureAtom<"SMARTS">(mol) -> Atom (null atom if there is no match)
    //

    template <ctll::fixed_string SMARTS>
    auto captureAtom(Molecule::Molecule auto &mol)
    {
        auto caps = capture<SMARTS>(mol);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // CTSmarts::capture<"SMARTS">(mol, atom) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS>
    auto capture(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        constexpr auto cap = captureMapping(smarts);
        auto [found, map] = iso.single(mol);
        return detail::captureMatchAtoms(mol, smarts, found, map, cap);
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, smarts, true, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol));
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique>
    auto captures(Molecule::Molecule auto &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, mapType};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol, atom) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, smarts, true, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol, atom));
    }

} // namespace Kitimar::CTSmarts
