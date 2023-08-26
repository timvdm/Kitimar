#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::match<"SMARTS">(mol) -> bool

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr bool match(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol)) // FIXME: use std::ranges::find_if -> check assembly?
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return true;
            return false;
        } else if constexpr (Config::specialize && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol))
                if (impl::singleBondMatch(smarts, mol, bond))
                    return true;
            return false;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, Config>{};
            return iso.match(mol);
        }
    }

    //
    // Atom
    //

    // ctse::match_atom<"SMARTS">(mol, atom) -> bool

    namespace impl {

        template<typename SmartsT, typename Config, Molecule::Molecule Mol>
        constexpr bool match_atom(Mol &mol, const auto &atom)
        {
            auto smarts = SmartsT{};
            if constexpr (smarts.isSingleAtom) {
                // Optimize single atom SMARTS
                return impl::singleAtomMatch(smarts, mol, atom);
            } else if constexpr (Config::specialize && smarts.isSingleBond) {
                // Optimize single bond SMARTS
                if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr))
                    return false;
                for (auto bond : get_bonds(mol, atom)) {
                    if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                        continue;
                    if (matchAtomExpr(mol, Kitimar::Molecule::get_nbr(mol, bond, atom), get<1>(smarts.atoms).expr))
                        return true;
                }
                return false;
            } else {
                auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, NoOptimizeConfig>{}; // FIXME: allow optimizations to be used...
                return iso.matchAtom(mol, atom);
            }
        }

    } // namespace impl

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr bool match_atom(Mol &mol, const auto &atom)
    {
        return impl::match_atom<Smarts<SMARTS>, Config, Mol>(mol, atom);

    }

    //
    // Bond
    //

    // ctse::match_bond<"SMARTS">(mol, bond) -> bool

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr bool match_bond(Mol &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr (Config::specialize && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            return impl::singleBondMatch(smarts, mol, bond);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::All, NoOptimizeConfig>{}; // FIXME: allow optimizations to be used...
            return iso.matchBond(mol, bond);
        }
    }

    //
    // Atom/Bond
    //

    // ctse::match<"SMARTS">(mol, atom) -> bool

    CTSMARTS_API_OVERLOAD_ATOM(match)

    // ctse::match<"SMARTS">(mol, bond) -> bool

    CTSMARTS_API_OVERLOAD_BOND(match)

} // mamespace Kitimar::CTSmarts
