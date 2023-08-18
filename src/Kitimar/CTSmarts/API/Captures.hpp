#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures(Mol &mol, MapTypeTag<M> mapType = {})
    {
        static_assert(M != MapType::Single, "Use CTSmarts::capture<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Smarts<SMARTS>{};
        static constexpr auto captureSet = captureMapping(smarts);
        //using Maps = IsomorphismMaps<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        using IndexMap = std::array<decltype(get_index(mol, get_atom(mol, 0))), captureSet.size() ? captureSet.size() : smarts.numAtoms>;
        using AtomMap = std::array<decltype(get_atom(mol, 0)), captureSet.size() ? captureSet.size() : smarts.numAtoms>;
        using AtomMaps = std::vector<AtomMap>;


        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            AtomMaps maps;
            maps.reserve(num_atoms(mol));
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    maps.push_back(impl::toCapture(mol, smarts, captureSet, true, IndexMap{get_index(mol, atom)}));
            return maps;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            AtomMaps maps;
            maps.reserve(num_bonds(mol));
            if constexpr (captureSet.size() == smarts.numAtoms) {
                for (auto bond : get_bonds(mol))
                    impl::singleBondCaptures(smarts, mol, bond, captureSet, maps, M == MapType::Unique);
            } else {
                for (auto bond : get_bonds(mol))
                    impl::singleBondCaptures(smarts, mol, bond, captureSet, maps, false);

                // Remove duplicates
                if constexpr (M == MapType::Unique) {
                    auto proj = [&mol] (const auto &atoms) { return impl::captureHash(mol, atoms); };
                    std::ranges::sort(maps, {}, proj);
                    const auto [first, last] = std::ranges::unique(maps, {}, proj);
                    maps.erase(first, last);
                }
            }
            return maps;
        } else {
            if constexpr (captureSet.size() == smarts.numAtoms) {
                auto iso = Isomorphism<Mol, decltype(smarts), M>{};
                if constexpr (__cpp_lib_ranges >= 202110L)
                    return iso.all(mol) | std::views::transform([&] (const auto &map) {
                        return impl::toCapture(mol, smarts, true, map, captureSet);
                    });
                else
                    // missing std::ranges::owning_view
                    return impl::toCaptures(mol, iso, captureSet, iso.all(mol));
            } else {
                auto iso = Isomorphism<Mol, decltype(smarts), MapType::All>{};
                AtomMaps maps = impl::toCaptures(mol, iso, captureSet, iso.all(mol));

                // Remove duplicates
                if constexpr (M == MapType::Unique) {
                    auto proj = [&mol] (const auto &atoms) { return impl::captureHash(mol, atoms); };
                    std::ranges::sort(maps, {}, proj);
                    const auto [first, last] = std::ranges::unique(maps, {}, proj);
                    maps.erase(first, last);
                }

                return maps;
            }
        }
    }

    // ctse::captures_unique<"SMARTS">(mol) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS,Molecule::Molecule Mol>
    constexpr auto captures_unique(Mol &mol)
    {
        return captures<SMARTS, MapType::Unique>(mol);
    }

    // ctse::captures_all<"SMARTS">(mol) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS,Molecule::Molecule Mol>
    constexpr auto captures_all(Mol &mol)
    {
        return captures<SMARTS, MapType::All>(mol);
    }

    //
    // Atom
    //

    // ctse::captures_atom<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures_atom(Mol &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(M != MapType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::capture_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        static constexpr auto captureSet = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.allAtom(mol, atom) | std::views::transform([&] (const auto &map) {
                return impl::toCapture(mol, smarts, captureSet, true, map);
            });
        else
            // missing std::ranges::owning_view
            return impl::toCaptures(mol, iso, captureSet, iso.allAtom(mol, atom));
    }

    // ctse::captures_atom_unique<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto captures_atom_unique(Molecule::Molecule auto &mol, const auto &atom)
    {
        return captures_atom<SMARTS, MapType::Unique>(mol, atom);
    }

    // ctse::captures_atom_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto captures_atom_all(Molecule::Molecule auto &mol, const auto &atom)
    {
        return captures_atom<SMARTS, MapType::All>(mol, atom);
    }

    //
    // Bond
    //

    // ctse::captures_bond<"SMARTS">(mol, bond, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures_bond(Mol &mol, const auto &bond, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(M != MapType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::capture_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        static constexpr auto captureSet = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.allBond(mol, bond) | std::views::transform([&] (const auto &map) {
                return impl::toCapture(mol, smarts, captureSet, true, map);
            });
        else
            // missing std::ranges::owning_view
            return impl::toCaptures(mol, iso, captureSet, iso.allBond(mol, bond));
    }

    // ctse::captures_bond_unique<"SMARTS">(mol, bond) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto captures_bond_unique(Molecule::Molecule auto &mol, const auto &bond)
    {
        return captures_bond<SMARTS, MapType::Unique>(mol, bond);
    }

    // ctse::captures_bond_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto captures_bond_all(Molecule::Molecule auto &mol, const auto &bond)
    {
        return captures_bond<SMARTS, MapType::All>(mol, bond);
    }


    //
    // Atom/Bond
    //

    // ctse::captures<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto captures(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, MapTypeTag<M> mapType = {})
    {
        return captures_atom<SMARTS>(mol, atom, mapType);
    }

    // ctse::captures<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto captures(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond, MapTypeTag<M> mapType = {})
    {
        return captures_bond<SMARTS>(mol, bond, mapType);
    }

    // ctse::captures_unique<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto captures_unique(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return captures_atom_unique<SMARTS>(mol, atom);
    }

    // ctse::captures_unique<"SMARTS">(mol, bond) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto captures_unique(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return captures_bond_unique<SMARTS>(mol, bond);
    }

    // ctse::captures_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto captures_all(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return captures_atom_all<SMARTS>(mol, atom);
    }

    // ctse::captures_all<"SMARTS">(mol, bond) -> std::vector<std::array<Atom, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto captures_all(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return captures_bond_all<SMARTS>(mol, bond);
    }


} // mamespace Kitimar::CTSmarts
