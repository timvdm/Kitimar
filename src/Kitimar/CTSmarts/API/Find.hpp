#pragma once

#include "Capture.hpp"

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::find<"SMARTS">(mol) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find(Molecule::Molecule auto &mol)
    {
        auto caps = capture<SMARTS, Config>(mol);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    // ctse::find_unique<"SMARTS">(mol) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find_unique(Molecule::Molecule auto &mol)
    {
        return captures_unique<SMARTS, Config>(mol) | std::views::transform([] (const auto &capture) {
            return capture[0];
        });
    }

    //
    // Atom
    //

    // ctse::find_atom<"SMARTS">(mol, atom) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find_atom(Molecule::Molecule auto &mol, const auto &atom)
    {
        auto caps = capture_atom<SMARTS, Config>(mol, atom);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // Bond
    //

    // ctse::find_bond<"SMARTS">(mol, bond) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find_bond(Molecule::Molecule auto &mol, const auto &bond)
    {
        auto caps = capture_bond<SMARTS, Config>(mol, bond);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // Atom/Bond
    //

    // ctse::find_atom<"SMARTS">(mol, atom) -> Atom (null atom if there is no match)

    CTSMARTS_API_OVERLOAD_ATOM(find)

    // ctse::find_bond<"SMARTS">(mol, bond) -> Atom (null atom if there is no match)

    CTSMARTS_API_OVERLOAD_BOND(find)

} // mamespace Kitimar::CTSmarts
