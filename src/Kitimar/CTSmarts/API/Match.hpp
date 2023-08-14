#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::match<"SMARTS">(mol) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr bool match(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol)) // FIXME: use std::ranges::find_if -> check assembly?
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return true;
            return false;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol))
                if (impl::singleBondMatch(smarts, mol, bond))
                    return true;
            return false;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single>{};
            return iso.match(mol);
        }
    }

    //
    // Atom
    //

    // ctse::match_atom<"SMARTS">(mol, atom) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr bool match_atom(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            return impl::singleAtomMatch(smarts, mol, atom);
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr))
                return false;
            for (auto bond : get_bonds(mol, atom)) {
                if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                    continue;
                if (matchAtomExpr(mol, get_nbr(mol, bond, atom), get<1>(smarts.atoms).expr))
                    return true;
            }
            return false;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single, NoOptimizationPolicy>{}; // FIXME: allow optimizations to be used...
            return iso.matchAtom(mol, atom);
        }
    }

    //
    // Bond
    //

    // ctse::match_bond<"SMARTS">(mol, bond) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr bool match_bond(Mol &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        //std::cout << "CTSmarts::bond<" << smarts.input() << ">(mol, " << get_index(mol, bond) << ")" << std::endl;
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            return impl::singleBondMatch(smarts, mol, bond);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::All, NoOptimizationPolicy>{}; // FIXME: allow optimizations to be used...
            return iso.matchBond(mol, bond);
        }
    }

    //
    // Atom/Bond
    //

    // ctse::match<"SMARTS">(mol, atom) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr bool match(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return match_atom<SMARTS>(mol, atom);
    }

    // ctse::match<"SMARTS">(mol, bond) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr bool match(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return match_bond<SMARTS>(mol, bond);
    }

} // mamespace Kitimar::CTSmarts
