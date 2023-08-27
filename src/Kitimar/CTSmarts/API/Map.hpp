#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::map<"SMARTS">(mol) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
        constexpr auto map(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        using Map = IsomorphismMap<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, Map{get_index(mol, atom)});
            return std::make_tuple(false, Map{});
        } else if constexpr (Config::specialize && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                auto m = impl::singleBondMap<Map>(smarts, mol, bond);
                if (std::get<0>(m))
                    return m;
            }
            return std::make_tuple(false, Map{});
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Molecule, Config>{};
            return iso.single(mol);
        }
    }

    //
    // Atom
    //

    // ctse::map_atom<"SMARTS">(mol, atom) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto map_atom(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Atom, Config>{};
        return iso.singleAtom(mol, atom);
    }

    //
    // Bond
    //

    // ctse::map_bond<"SMARTS">(mol, bond) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto map_bond(Mol &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Bond, Config>{};
        return iso.singleBond(mol, bond);
    }

    //
    // Atom/Bond
    //

    // ctse::map<"SMARTS">(mol, atom) -> std::tuple<bool, std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM(map)

    // ctse::map<"SMARTS">(mol, bond) -> std::tuple<bool, std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND(map)

} // mamespace Kitimar::CTSmarts
