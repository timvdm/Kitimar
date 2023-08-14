#pragma once

#include "API/Match.hpp"
#include "API/Count.hpp"
#include "API/Map.hpp"
#include "API/Maps.hpp"


namespace Kitimar::CTSmarts {



    /*
     *
     * bool match(mol)
     * bool match_atom(mol, atom)
     * bool match_bond(mol, bond)
     * bool match(mol, atom/bond)
     *
     *
     * int count(mol, type = Unique)
     * int count_unique(mol)
     * int count_all(mol)
     *
     * int count_atom(mol, atom, type = Unique)
     * int count_atom_unique(mol, atom)
     * int count_atom_all(mol, atom)
     *
     * int count_bond(mol, bond, type = Unique)
     * int count_bond_unique_bond(mol, bond)
     * int count_bond_all(mol, bond)
     *
     * int count(mol, atom/bond, type = Unique)
     * int count_unique(mol, atom/bond)
     * int count_all(mol, atom/bond)
     *
     *
     * Map map(mol)
     * Map map_atom(mol, atom)
     * Map map_bond(mol, bond)
     * Map map(mol, atom/bond)
     *
     *
     *
     *
     *
     *
     *
     * Maps maps(mol, type = Unique)
     * Maps maps_unique(mol)
     * Maps maps_all(mol)
     *
     * Maps maps_atom(mol, atom, type = Unqiue)
     * Maps maps_atom_unique(mol, atom)
     * Maps maps_atom_all(mol, atom)
     *
     * Maps maps_bond(mol, bond, type = Unique)
     * Maps maps_bond_unique(mol, bond)
     * Maps maps_bond_all(mol, bond)
     *
     * Maps maps(mol, atom/bond, type = Unique)
     * Maps maps_unique(mol, atom/bond)
     * Maps maps_all(mol, atom/bond)
     *
     *
     *
     * (bool, Atom...) capture(mol)
     * (bool, Atom...) capture(mol, atom/bond)
     *     (bool, Atom...) capture_atom(mol, atom)
     *     (bool, Atom...) capture_bond(mol, bond)
     *
     * Maps captures(mol, type = Unique)
     * Maps captures(mol, atom/bond, type = Unique)
     *     Maps captures_atom(mol, atom, type = Unique)
     *     Maps captures_bond(mol, bond, type = Unique)
     *
     * Maps captures_unique(mol)
     * Maps captures_unique(mol, atom/bond)
     *     Maps captures_atom_unique(mol, atom)
     *     Maps captures_bond_unique(mol)
     *
     * Maps captures_all(mol)
     * Maps captures_all(mol, atom/bond)
     *     Maps captures_atom_all(mol, atom)
     *     Maps captures_bond_all(mol)
     *
     *
     *
     */









    //
    // CTSmarts::capture<"SMARTS">(mol) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    auto capture(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, atom);
            return std::make_tuple(false, null_atom(mol));
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            constexpr auto cap = captureMapping(smarts);
            for (auto bond : get_bonds(mol)) {
                auto matchType = impl::singleBondMatch(smarts, mol, bond);
                if (matchType)
                    return impl::singleBondCapture(smarts, mol, bond, matchType);
            }
            return impl::singleBondCapture(smarts, mol, null_bond(mol), 0);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single>{};
            constexpr auto cap = captureMapping(smarts);
            auto [found, map] = iso.single(mol);
            return impl::captureMatchAtoms(mol, smarts, found, map, cap);
        }
    }

    //
    // CTSmarts::captureAtom<"SMARTS">(mol) -> Atom (null atom if there is no match)
    //

    template <ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    auto captureAtom(Mol &mol)
    {
        auto caps = capture<SMARTS>(mol);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // CTSmarts::capture<"SMARTS">(mol, atom) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    auto capture(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        constexpr auto cap = captureMapping(smarts);
        auto [found, map] = iso.single(mol);
        return impl::captureMatchAtoms(mol, smarts, found, map, cap);
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures(Mol &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism<Mol, decltype(smarts), M>{};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol) | std::views::transform([&] (const auto &map) {
                return impl::captureAtoms(mol, smarts, true, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return impl::copyCapture(mol, iso, cap, iso.all(mol));
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures(Mol &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism<Mol, decltype(smarts), M>{};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol, atom) | std::views::transform([&] (const auto &map) {
                return impl::captureAtoms(mol, smarts, true, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return impl::copyCapture(mol, iso, cap, iso.all(mol, atom));
    }

} // namespace Kitimar::CTSmarts

namespace ctse = Kitimar::CTSmarts;
