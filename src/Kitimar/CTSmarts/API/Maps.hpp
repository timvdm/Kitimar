#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::maps<"SMARTS">(mol, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto maps(const Mol &mol, SearchTypeTag<M> = {})
    {
        static_assert(M != SearchType::Single, "Use CTSmarts::map<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
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
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            Maps maps;
            maps.reserve(num_bonds(mol));
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == SearchType::Unique) {
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
            auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Molecule, Config>{};
            return iso.all(mol);
        }
    }

    // ctse::maps_unique<"SMARTS">(mol) -> std::vector<std::array<int, N>>

    CTSMARTS_API_UNIQUE(maps)

    // ctse::maps_all<"SMARTS">(mol) -> std::vector<std::array<int, N>>

    CTSMARTS_API_ALL(maps)

    //
    // Atom
    //

    // ctse::maps_atom<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto maps_atom(const Mol &mol, const auto &atom, SearchTypeTag<M> SearchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::map_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Atom, Config>{};
        return iso.allAtom(mol, atom);
    }

    // ctse::maps_atom_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_ATOM_UNIQUE(maps)

    // ctse::maps_atom_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_ATOM_ALL(maps)

    //
    // Bond
    //

    // ctse::maps_bond<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto maps_bond(const Mol &mol, const auto &bond, SearchTypeTag<M> SearchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom /*&& !smarts.isSingleBond*/,
                "Use CTSmarts::map_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Bond, Config>{};
        return iso.allBond(mol, bond);
    }

    // ctse::maps_bond_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_BOND_UNIQUE(maps)

    // ctse::maps_bond_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_BOND_ALL(maps)

    //
    // Atom/Bond
    //

    // ctse::maps<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_SEARCH(maps)

    // ctse::maps<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_SEARCH(maps)

    // ctse::maps_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(maps)

    // ctse::maps_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_UNIQUE(maps)

    // ctse::maps_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_ALL(maps)

    // ctse::maps_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_ALL(maps)

} // mamespace Kitimar::CTSmarts
