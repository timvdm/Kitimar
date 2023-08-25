#pragma once

#include "Match.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::count<"SMARTS">(mol, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto count(Mol &mol, SearchTypeTag<M> = {})
    {
        static_assert(M != SearchType::Single, "Use CTSmarts::contains<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            auto n = 0;
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    ++n;
            return n;
        } else if constexpr (Config::specialize && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            auto n = 0;
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == SearchType::Unique) {
                    if (impl::singleBondMatch(smarts, mol, bond))
                        ++n;
                } else {
                    n += impl::singleBondCount(smarts, mol, bond);
                }
            }
            return n;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), M, Config>{};
            return iso.count(mol);
        }
    }

    // ctse::count_unique<"SMARTS">(mol) -> std::integeral

    CTSMARTS_API_UNIQUE(count)

    // ctse::count_all<"SMARTS">(mol) -> std::integeral

    CTSMARTS_API_ALL(count)

    //
    // Atom
    //

    // ctse::count_atom<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto count_atom(Mol &mol, const auto &atom, SearchTypeTag<M> = {})
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom)
            return match_atom<SMARTS, Config>(mol, atom) ? 1 : 0;
        else {
            auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizeConfig>{}; // FIXME: NoOptimizeConfig
            return iso.countAtom(mol, atom);
        }
    }

    // ctse::count_atom_unique<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_ATOM_UNIQUE(count)

    // ctse::count_atom_all<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_ATOM_ALL(count)

    //
    // Bond
    //

    // ctse::count_bond<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto count_bond(Mol &mol, const auto &bond, SearchTypeTag<M> = {})
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr (Config::specialize && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            if constexpr (M == SearchType::Unique) {
                if (impl::singleBondMatch(smarts, mol, bond))
                    return 1;
            } else {
                return impl::singleBondCount(smarts, mol, bond);
            }
            return 0;
        }
        auto iso = Isomorphism<Mol, decltype(smarts), M, NoOptimizeConfig>{}; // FIXME: NoOptimizeConfig
        return iso.countBond(mol, bond);
    }

    // ctse::count_bond_unique<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_BOND_UNIQUE(count)

    // ctse::count_bond_all<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_BOND_ALL(count)

    //
    // Atom/Bond
    //

    // ctse::count<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::integeral

    CTSMARTS_API_OVERLOAD_ATOM_SEARCH(count)

    // ctse::count<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::integeral

    CTSMARTS_API_OVERLOAD_BOND_SEARCH(count)

    // ctse::count_unique<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(count)

    // ctse::count_unique<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_OVERLOAD_BOND_UNIQUE(count)

    // ctse::count_all<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_OVERLOAD_ATOM_ALL(count)

    // ctse::count_all<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_OVERLOAD_BOND_ALL(count)

} // mamespace Kitimar::CTSmarts
