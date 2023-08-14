#pragma once

#include "Isomorphism.hpp"

#include <Kitimar/Molecule/Molecule.hpp>

namespace Kitimar::CTSmarts {

    namespace detail {

        template<typename Smarts, auto NumCaptures>
        auto captureAtoms(Molecule::Molecule auto &mol, Smarts, bool found, const auto &map,
                          const std::array<int, NumCaptures> &cap)
        {
            using Atom = decltype(get_atom(mol, 0));
            if constexpr (NumCaptures) {
                std::array<Atom, NumCaptures> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < NumCaptures; ++i)
                        atoms[i] = get_atom(mol, map[cap[i]]);
                return atoms;
            } else {
                std::array<Atom, Smarts::numAtoms> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0; i < Smarts::numAtoms; ++i)
                        atoms[i] = get_atom(mol, map[i]); // FIXME: null atoms...
                return atoms;
            }
        }

        auto captureMatchAtoms(Molecule::Molecule auto &mol, auto smarts, bool found, const auto &map, const auto &cap)
        {
            return std::tuple_cat(std::make_tuple(found), captureAtoms(mol, smarts, found, map, cap));
        }


        template<auto N>
        auto copyCapture(Molecule::Molecule auto &mol, const auto &iso, const std::array<int, N> &cap, const auto &caps) noexcept
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto M = N ? N : iso.smarts.numAtoms;
            std::vector<std::array<Atom, M>> v;
            auto r = caps | std::views::transform([&] (const auto &map) {
                return captureAtoms(mol, iso.smarts, true, map, cap);
            });
            std::ranges::copy(r, std::back_inserter(v));
            return v;
        }

        constexpr bool singleAtomMatch(auto smarts, auto &mol, const auto &atom)
        {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr);
        }

        constexpr bool singleBondMatchHelper(auto smarts, auto &mol, const auto &bond, const auto &source, const auto &target)
        {
            return matchAtomExpr(mol, source, get<0>(smarts.atoms).expr) && matchAtomExpr(mol, target, get<1>(smarts.atoms).expr);
        }

        // 0 -> no match
        // 1 -> source is SMARTS atom 0, target is SMARTS atom 1
        // 2 -> source is SMARTS atom 1, target is SMARTS atom 0
        constexpr int singleBondMatch(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return 1;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return 2;
            return 0;
        }

        constexpr int singleBondCount(auto smarts, auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            auto n = 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                ++n;;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                ++n;
            return n;
        }

        constexpr auto singleBondCapture(auto smarts, auto &mol, const auto &bond, int singleBondMatchType)
        {
            constexpr auto cap = captureMapping(smarts);
            if constexpr (cap.size() == 1) {
                switch (singleBondMatchType) {
                    case 1:
                        return cap[0] == 0 ? std::make_tuple(true, get_source(mol, bond)) :
                                             std::make_tuple(true, get_target(mol, bond));
                    case 2:
                        return cap[0] == 0 ? std::make_tuple(true, get_target(mol, bond)) :
                                             std::make_tuple(true, get_source(mol, bond));
                    default:
                        return std::make_tuple(false, null_atom(mol));
                }
            } else {
                if (!singleBondMatchType)
                    return std::make_tuple(false, null_atom(mol), null_atom(mol));

                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                if constexpr (cap.size())
                    if (cap[0] > cap[1])
                        std::swap(source, target);

                if (singleBondMatchType == 1)
                    return std::make_tuple(true, source, target);
                return std::make_tuple(true, target, source);
            }

        }

        constexpr auto centralAtomMap(auto smarts, auto &mol, const auto &atom)
        {
            std::array<int, smarts.numAtoms> map;
        }


    } // namespace detail

    /*
     *
     * bool match(mol)
     * bool match_atom(mol, atom)
     * bool match_bond(mol, bond)
     * bool match(mol, atom/bond)
     *
     *
     * int count(mol, type = Unique)
     * int count_unique(mol)
     * int count_all(mol)
     *
     * int count_atom(mol, atom, type = Unique)
     * int count_atom_unique(mol, atom)
     * int count_atom_all(mol, atom)
     *
     * int count_bond(mol, bond, type = Unique)
     * int count_bond_unique_bond(mol, bond)
     * int count_bond_all(mol, bond)
     *
     * int count(mol, atom/bond, type = Unique)
     * int count_unique(mol, atom/bond)
     * int count_all(mol, atom/bond)
     *
     *
     * Map map(mol)
     * Map map_atom(mol, atom)
     * Map map_bond(mol, bond)
     * Map map(mol, atom/bond)
     *
     *
     *
     *
     *
     *
     *
     * Maps maps(mol, type = Unique)
     * Maps maps_unique(mol)
     * Maps maps_all(mol)
     *
     * Maps maps_atom(mol, atom, type = Unqiue)
     * Maps maps_atom_unique(mol, atom)
     * Maps maps_atom_all(mol, atom)
     *
     * Maps maps_bond(mol, bond, type = Unique)
     * Maps maps_bond_unique(mol, bond)
     * Maps maps_bond_all(mol, bond)
     *
     * Maps maps(mol, atom/bond, type = Unique)
     * Maps maps_unique(mol, atom/bond)
     * Maps maps_all(mol, atom/bond)
     *
     *
     *
     * (bool, Atom...) capture(mol)
     * (bool, Atom...) capture(mol, atom/bond)
     *     (bool, Atom...) capture_atom(mol, atom)
     *     (bool, Atom...) capture_bond(mol, bond)
     *
     * Maps captures(mol, type = Unique)
     * Maps captures(mol, atom/bond, type = Unique)
     *     Maps captures_atom(mol, atom, type = Unique)
     *     Maps captures_bond(mol, bond, type = Unique)
     *
     * Maps captures_unique(mol)
     * Maps captures_unique(mol, atom/bond)
     *     Maps captures_atom_unique(mol, atom)
     *     Maps captures_bond_unique(mol)
     *
     * Maps captures_all(mol)
     * Maps captures_all(mol, atom/bond)
     *     Maps captures_atom_all(mol, atom)
     *     Maps captures_bond_all(mol)
     *
     *
     *
     */



    //
    // Match
    //

    // ctse::match<"SMARTS">(mol) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr bool match(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol)) // FIXME: use std::ranges::find_if -> check assembly?
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return true;
            return false;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol))
                if (detail::singleBondMatch(smarts, mol, bond))
                    return true;
            return false;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single>{};
            return iso.match(mol);
        }
    }

    // ctse::match_atom<"SMARTS">(mol, atom) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr bool match_atom(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            return detail::singleAtomMatch(smarts, mol, atom);
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

    // ctse::match_bond<"SMARTS">(mol, bond) -> bool

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr bool match_bond(Mol &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        //std::cout << "CTSmarts::bond<" << smarts.input() << ">(mol, " << get_index(mol, bond) << ")" << std::endl;
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            return detail::singleBondMatch(smarts, mol, bond);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::All, NoOptimizationPolicy>{}; // FIXME: allow optimizations to be used...
            return iso.matchBond(mol, bond);
        }
    }

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

    //
    // Count
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
                if (detail::singleAtomMatch(smarts, mol, atom))
                    ++n;
            return n;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            auto n = 0;
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == MapType::Unique) {
                    if (detail::singleBondMatch(smarts, mol, bond))
                        ++n;
                } else {
                    n += detail::singleBondCount(smarts, mol, bond);
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

    //
    // Map
    //

    // ctse::map<"SMARTS">(mol) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
        constexpr auto map(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        using Map = IsomorphismMap<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, Map{get_index(mol, atom)});
            return std::make_tuple(false, Map{});
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                switch (detail::singleBondMatch(smarts, mol, bond)) {
                    case 1:
                        return std::make_tuple(true, Map{get_index(mol, get_source(mol, bond)),
                                                         get_index(mol, get_target(mol, bond))});
                    case 2:
                        return std::make_tuple(true, Map{get_index(mol, get_target(mol, bond)),
                                                         get_index(mol, get_source(mol, bond))});
                    default:
                        break;
                }
            }
            return std::make_tuple(false, Map{});
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single>{};
            return iso.single(mol);
        }
    }

    // ctse::map_atom<"SMARTS">(mol, atom) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr auto map_atom(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        return iso.singleAtom(mol, atom);
    }

    // ctse::map<"SMARTS">(mol, atom) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
        requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto map(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Atom atom)
    {
        return map_atom<SMARTS>(mol, atom);
    }

    // ctse::map_bond<"SMARTS">(mol, bond) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    constexpr auto map_bond(Mol &mol, const auto &bond)
    {
        auto smarts = Smarts<SMARTS>{};
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single, NoOptimizationPolicy>{}; // FIXME: NoOptimizationPolicy
        return iso.singleBond(mol, bond);
    }

    // ctse::map<"SMARTS">(mol, bond) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, Molecule::Molecule Mol>
        requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType)
    constexpr auto map(Mol &mol, typename Molecule::MoleculeTraits<Mol>::Bond bond)
    {
        return map_bond<SMARTS>(mol, bond);
    }

    //
    // Maps
    //

    //
    // CTSmarts::multi<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::vector<std::vector<int>>
    //

    template<ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    constexpr auto multi(Mol &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        using Maps = IsomorphismMaps<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        using Map = Maps::value_type;

        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            Maps maps;
            maps.reserve(num_atoms(mol));
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    maps.push_back(Map{get_index(mol, atom)});
            return maps;
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            Maps maps;
            maps.reserve(num_bonds(mol));
            for (auto bond : get_bonds(mol)) {
                switch (detail::singleBondMatch(smarts, mol, bond)) {
                    case 1:
                        maps.push_back(Map{get_index(mol, get_source(mol, bond)), get_index(mol, get_target(mol, bond))});
                        if constexpr (M == MapType::All)
                            maps.push_back(Map{get_index(mol, get_target(mol, bond)), get_index(mol, get_source(mol, bond))});
                        break;
                    case 2:
                        maps.push_back(Map{get_index(mol, get_target(mol, bond)), get_index(mol, get_source(mol, bond))});
                        if constexpr (M == MapType::All)
                            maps.push_back(Map{get_index(mol, get_source(mol, bond)), get_index(mol, get_target(mol, bond))});

                        break;
                    default:
                        break;
                }
            }
            return maps;
        //} else if constexpr (smarts.centralAtom != -1) {

        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), M>{};
            return iso.all(mol);
        }
    }

    //
    // CTSmarts::capture<"SMARTS">(mol) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    auto capture(Mol &mol)
    {
        auto smarts = Smarts<SMARTS>{};
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (detail::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, atom);
            return std::make_tuple(false, null_atom(mol));
        } else if constexpr (smarts.isSingleBond) {
            // Optimize single bond SMARTS
            constexpr auto cap = captureMapping(smarts);
            for (auto bond : get_bonds(mol)) {
                auto matchType = detail::singleBondMatch(smarts, mol, bond);
                if (matchType)
                    return detail::singleBondCapture(smarts, mol, bond, matchType);
            }
            return detail::singleBondCapture(smarts, mol, null_bond(mol), 0);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), MapType::Single>{};
            constexpr auto cap = captureMapping(smarts);
            auto [found, map] = iso.single(mol);
            return detail::captureMatchAtoms(mol, smarts, found, map, cap);
        }
    }

    //
    // CTSmarts::captureAtom<"SMARTS">(mol) -> Atom (null atom if there is no match)
    //

    template <ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    auto captureAtom(Mol &mol)
    {
        auto caps = capture<SMARTS>(mol);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // CTSmarts::capture<"SMARTS">(mol, atom) -> tuple<bool, Atom...>
    //

    template <ctll::fixed_string SMARTS, Molecule::Molecule Mol>
    auto capture(Mol &mol, const auto &atom)
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism{smarts, Single};
        constexpr auto cap = captureMapping(smarts);
        auto [found, map] = iso.single(mol);
        return detail::captureMatchAtoms(mol, smarts, found, map, cap);
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures(Mol &mol, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism<Mol, decltype(smarts), M>{};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, smarts, true, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol));
    }

    //
    // CTSmarts::captures<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::range<std::array<Atom, N>>
    //

    template <ctll::fixed_string SMARTS, MapType M = MapType::Unique, Molecule::Molecule Mol>
    auto captures(Mol &mol, const auto &atom, MapTypeTag<M> mapType = {})
    {
        auto smarts = Smarts<SMARTS>{};
        auto iso = Isomorphism<Mol, decltype(smarts), M>{};
        static constexpr auto cap = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.all(mol, atom) | std::views::transform([&] (const auto &map) {
                return detail::captureAtoms(mol, smarts, true, map, cap);
            });
        else
            // missing std::ranges::owning_view
            return detail::copyCapture(mol, iso, cap, iso.all(mol, atom));
    }

} // namespace Kitimar::CTSmarts

namespace ctse = Kitimar::CTSmarts;
