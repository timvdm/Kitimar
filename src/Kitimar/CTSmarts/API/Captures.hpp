#pragma once

#include "Util.hpp"
#include "../Isomorphism.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto captures(const Mol &mol, SearchTypeTag<M> = {})
    {
        static_assert(M != SearchType::Single, "Use CTSmarts::capture<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static constexpr auto captureSet = captureMapping(smarts);
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
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            AtomMaps maps;
            maps.reserve(num_bonds(mol));
            if constexpr (captureSet.size() == smarts.numAtoms) {
                for (auto bond : get_bonds(mol))
                    impl::singleBondCaptures(smarts, mol, bond, captureSet, maps, M == SearchType::Unique);
            } else {
                for (auto bond : get_bonds(mol))
                    impl::singleBondCaptures(smarts, mol, bond, captureSet, maps, false);

                // Remove duplicates
                if constexpr (M == SearchType::Unique) {
                    auto proj = [&mol] (const auto &atoms) { return impl::captureHash(mol, atoms); };
                    std::ranges::sort(maps, {}, proj);
                    const auto [first, last] = std::ranges::unique(maps, {}, proj);
                    maps.erase(first, last);
                }
            }
            return maps;
        } else {
            if constexpr (captureSet.size() == smarts.numAtoms) {
                auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Molecule, Config>{};
                if constexpr (__cpp_lib_ranges >= 202110L)
                    return iso.all(mol) | std::views::transform([&] (const auto &map) {
                        return impl::toCapture(mol, smarts, true, map, captureSet);
                    });
                else
                    // missing std::ranges::owning_view
                    return impl::toCaptures(mol, iso, captureSet, iso.all(mol));
            } else {
                auto iso = Isomorphism<Mol, decltype(smarts), SearchType::All, SeedType::Molecule, Config>{};
                AtomMaps maps = impl::toCaptures(mol, iso, captureSet, iso.all(mol));

                // Remove duplicates
                if constexpr (M == SearchType::Unique) {
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

    CTSMARTS_API_UNIQUE(captures)

    // ctse::captures_all<"SMARTS">(mol) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_ALL(captures)

    //
    // Atom
    //

    // ctse::captures_atom<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto captures_atom(const Mol &mol, const auto &atom, SearchTypeTag<M> searchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::capture_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Atom, Config>{};
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

    CTSMARTS_API_ATOM_UNIQUE(captures)

    // ctse::captures_atom_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_ATOM_ALL(captures)

    //
    // Bond
    //

    // ctse::captures_bond<"SMARTS">(mol, bond, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, SearchType M = SearchType::Unique, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto captures_bond(const Mol &mol, const auto &bond, SearchTypeTag<M> searchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::capture_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Bond, Config>{};
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

    CTSMARTS_API_BOND_UNIQUE(captures)

    // ctse::captures_bond_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_BOND_ALL(captures)

    //
    // Atom/Bond
    //

    // ctse::captures<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_SEARCH(captures)

    // ctse::captures<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_SEARCH(captures)

    // ctse::captures_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(captures)

    // ctse::captures_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_UNIQUE(captures)

    // ctse::captures_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_ALL(captures)

    // ctse::captures_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_ALL(captures)

} // mamespace Kitimar::CTSmarts
