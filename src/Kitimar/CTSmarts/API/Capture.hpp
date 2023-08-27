#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::capture<"SMARTS">(mol) -> std::tuple<bool, Atom...>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto capture(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, atom);
            return std::make_tuple(false, null_atom(mol));
        } else if constexpr (Config::specialize && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                auto matchType = impl::singleBondMatch(smarts, mol, bond);
                if (matchType)
                    return impl::singleBondCapture(smarts, mol, bond, matchType);
            }
            return impl::singleBondCapture(smarts, mol, null_bond(mol), 0);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Molecule, Config>{};
            constexpr auto captureSet = captureMapping(smarts);
            auto [found, map] = iso.single(mol);
            return impl::captureMatchAtoms(mol, smarts, captureSet, found, map);
        }
    }

    //
    // Atom
    //

    // ctse::capture_atom<"SMARTS">(mol, atom) -> std::tuple<bool, Atom...>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto capture_atom(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Atom, Config>{};
        constexpr auto captureSet = captureMapping(smarts);
        auto [found, map] = iso.singleAtom(mol, atom);
        return impl::captureMatchAtoms(mol, smarts, captureSet, found, map);
    }

    //
    // Bond
    //

    // ctse::capture_bond<"SMARTS">(mol, bond) -> std::tuple<bool, Atom...>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto capture_bond(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Bond, Config>{};
        constexpr auto captureSet = captureMapping(smarts);
        auto [found, map] = iso.singleBond(mol, atom);
        return impl::captureMatchAtoms(mol, smarts, captureSet, found, map);
    }

    //
    // Atom/Bond
    //

    // ctse::capture<"SMARTS">(mol, atom) -> std::tuple<bool, Atom...>

    CTSMARTS_API_OVERLOAD_ATOM(capture)

    // ctse::capture<"SMARTS">(mol, bond) -> std::tuple<bool, Atom...>

    CTSMARTS_API_OVERLOAD_BOND(capture)


} // mamespace Kitimar::CTSmarts
