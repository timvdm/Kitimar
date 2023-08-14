#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::maps<"SMARTS">(mol, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto maps(Mol &mol, MapTypeTag<M> mapType = {})
    {
        static_assert(M != MapType::Single, "Use CTSmarts::map<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Smarts<SMARTS>{};
        using Maps = IsomorphismMaps<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        using Map = Maps::value_type;

        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            Maps maps;
            maps.reserve(num_atoms(mol));
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    maps.push_back(Map{get_index(mol, atom)});
            return maps;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            Maps maps;
            maps.reserve(num_bonds(mol));
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == MapType::Unique) {
                    auto [found, map] = impl::singleBondMap<Map>(smarts, mol, bond);
                    if (found)
                        maps.push_back(map);
                } else {
                    impl::singleBondMaps(smarts, mol, bond, maps);
                }
            }
            return maps;
        //} else if constexpr (smarts.centralAtom != -1) {

        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), M>{};
            return iso.all(mol);
        }
    }

    // ctse::maps_unique<"SMARTS">(mol) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS,Molecule::Molecule Mol>
    constexpr auto maps_unique(Mol &mol)
    {
        return maps<SMARTS, MapType::Unique>(mol);
    }

    // ctse::maps_all<"SMARTS">(mol) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS,Molecule::Molecule Mol>
    constexpr auto maps_all(Mol &mol)
    {
        return maps<SMARTS, MapType::All>(mol);
    }

    //
    // Atom
    //

    // ctse::maps_atom<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto maps_atom(Mol &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(M != MapType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::map_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        return iso.allAtom(mol, atom);
    }

    // ctse::maps_atom_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto maps_atom_unique(Molecule::Molecule auto &mol, const auto &atom)
    {
        return maps_atom<SMARTS, MapType::Unique>(mol, atom);
    }

    // ctse::maps_atom_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto maps_atom_all(Molecule::Molecule auto &mol, const auto &atom)
    {
        return maps_atom<SMARTS, MapType::All>(mol, atom);
    }

    //
    // Bond
    //

    // ctse::maps_bond<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto maps_bond(Mol &mol, const auto &bond, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(M != MapType::Single && !smarts.isSingleAtom && !smarts.isSingleBond,
                "Use CTSmarts::map_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        return iso.allBond(mol, bond);
    }

    // ctse::maps_bond_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto maps_bond_unique(Molecule::Molecule auto &mol, const auto &bond)
    {
        return maps_bond<SMARTS, MapType::Unique>(mol, bond);
    }

    // ctse::maps_bond_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS>
    constexpr auto maps_bond_all(Molecule::Molecule auto &mol, const auto &bond)
    {
        return maps_bond<SMARTS, MapType::All>(mol, bond);
    }

    //
    // Atom/Bond
    //

    // ctse::maps<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto maps(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, MapTypeTag<M> mapType = {})
    {
        return maps_atom<SMARTS>(mol, atom, mapType);
    }

    // ctse::maps<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto maps(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond, MapTypeTag<M> mapType = {})
    {
        return maps_bond<SMARTS>(mol, bond, mapType);
    }

    // ctse::maps_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto maps_unique(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return maps_atom_unique<SMARTS>(mol, atom);
    }

    // ctse::maps_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto maps_unique(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return maps_bond_unique<SMARTS>(mol, bond);
    }

    // ctse::maps_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto maps_all(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return maps_atom_all<SMARTS>(mol, atom);
    }

    // ctse::maps_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto maps_all(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return maps_bond_all<SMARTS>(mol, bond);
    }

} // mamespace Kitimar::CTSmarts
