#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::count<"SMARTS">(mol, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto count(Mol &mol, MapTypeTag<M> mapType = {})
    {
        static_assert(M != MapType::Single, "Use CTSmarts::contains<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            auto n = 0;
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    ++n;
            return n;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            auto n = 0;
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == MapType::Unique) {
                    if (impl::singleBondMatch(smarts, mol, bond))
                        ++n;
                } else {
                    n += impl::singleBondCount(smarts, mol, bond);
                }
            }
            return n;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), M>{};
            return iso.count(mol);
        }
    }

    // ctse::count_unique<"SMARTS">(mol) -> std::integeral

    template<ctll::fixed_string SMARTS>
    constexpr auto count_unique(Molecule::Molecule auto &mol)
    {
        return count<SMARTS, MapType::Unique>(mol);
    }

    // ctse::count_all<"SMARTS">(mol) -> std::integeral

    template<ctll::fixed_string SMARTS>
    constexpr auto count_all(Molecule::Molecule auto &mol)
    {
        return count<SMARTS, MapType::All>(mol);
    }

    //
    // Atom
    //

    // ctse::count_atom<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto count_atom(Mol &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(M != MapType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        return iso.countAtom(mol, atom);
    }

    // ctse::count_atom_unique<"SMARTS">(mol, atom) -> std::integeral

    template<ctll::fixed_string SMARTS>
    constexpr auto count_atom_unique(Molecule::Molecule auto &mol, const auto &atom)
    {
        return count_atom<SMARTS, MapType::Unique>(mol, atom);
    }

    // ctse::count_atom_all<"SMARTS">(mol, atom) -> std::integeral

    template<ctll::fixed_string SMARTS>
    constexpr auto count_atom_all(Molecule::Molecule auto &mol, const auto &atom)
    {
        return count_atom<SMARTS, MapType::All>(mol, atom);
    }

    //
    // Bond
    //

    // ctse::count_bond<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto count_bond(Mol &mol, const auto &bond, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(M != MapType::Single && !smarts.isSingleAtom && !smarts.isSingleBond,
                "Use CTSmarts::match_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        return iso.countBond(mol, bond);
    }

    // ctse::count_bond_unique<"SMARTS">(mol, bond) -> std::integeral

    template<ctll::fixed_string SMARTS>
    constexpr auto count_bond_unique(Molecule::Molecule auto &mol, const auto &bond)
    {
        return count_bond<SMARTS, MapType::Unique>(mol, bond);
    }

    // ctse::count_bond_all<"SMARTS">(mol, bond) -> std::integeral

    template<ctll::fixed_string SMARTS>
    constexpr auto count_bond_all(Molecule::Molecule auto &mol, const auto &bond)
    {
        return count_bond<SMARTS, MapType::All>(mol, bond);
    }

    //
    // Atom/Bond
    //

    // ctse::count<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto count(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom, MapTypeTag<M> mapType = {})
    {
        return count_atom<SMARTS>(mol, atom, mapType);
    }

    // ctse::count<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto count(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond, MapTypeTag<M> mapType = {})
    {
        return count_bond<SMARTS>(mol, bond, mapType);
    }

    // ctse::count_unique<"SMARTS">(mol, atom) -> std::integeral

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto count_unique(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return count_atom_unique<SMARTS>(mol, atom);
    }

    // ctse::count_unique<"SMARTS">(mol, bond) -> std::integeral

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto count_unique(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return count_bond_unique<SMARTS>(mol, bond);
    }

    // ctse::count_all<"SMARTS">(mol, atom) -> std::integeral

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto count_all(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return count_atom_all<SMARTS>(mol, atom);
    }

    // ctse::count_all<"SMARTS">(mol, bond) -> std::integeral

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto count_all(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return count_bond_all<SMARTS>(mol, bond);
    }

} // mamespace Kitimar::CTSmarts
