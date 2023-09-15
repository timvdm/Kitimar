#pragma once

#include "API/Match.hpp"
#include "API/Count.hpp"
#include "API/Map.hpp"
#include "API/Maps.hpp"
#include "API/Capture.hpp"
#include "API/Captures.hpp"
#include "API/Find.hpp"

/*

API
===

Match
-----

    match(mol)
    match_atom(mol, atom)
    match_bond(mol, bond)
    match(mol, atom/bond)

    Return type: bool

Count
-----

    count(mol, type = Unique)
    count_unique(mol)
    count_all(mol)

    count_atom(mol, atom, type = Unique)
    count_atom_unique(mol, atom)
    count_atom_all(mol, atom)

    count_bond(mol, bond, type = Unique)
    count_bond_unique_bond(mol, bond)
    count_bond_all(mol, bond)

    count(mol, atom/bond, type = Unique)
    count_unique(mol, atom/bond)
    count_all(mol, atom/bond)

    Return type: int

Map
---

    map(mol)
    map_atom(mol, atom)
    map_bond(mol, bond)
    map(mol, atom/bond)

    Return type: std::tuple<bool, std::array<int, N>>

    Example
    ~~~~~~~

        auto [found, map] = ctse::map<"C=O">(mol);
        if (found) {
            // Use map...
        }

Maps
----

    maps(mol, type = Unique)
    maps_unique(mol)
    maps_all(mol)

    maps_atom(mol, atom, type = Unqiue)
    maps_atom_unique(mol, atom)
    maps_atom_all(mol, atom)

    maps_bond(mol, bond, type = Unique)
    maps_bond_unique(mol, bond)
    maps_bond_all(mol, bond)

    maps(mol, atom/bond, type = Unique)
    maps_unique(mol, atom/bond)
    maps_all(mol, atom/bond)

    Return type: std::vector<std::array<int, N>>

Capture
-------

    capture(mol)
    capture_atom(mol, atom)
    capture_bond(mol, bond)
    capture(mol, atom/bond)

    Return type: std::tuple<bool, Atom...>

    Examples
    ~~~~~~~~

        auto [found, C, O] = ctse::capture<"C=O">
        if (found) {
            // Use atoms C and O...
        }

        auto [found, C, O] = ctse::capture<"N[C:1]=[O:2]">
        if (found) {
            // Use atoms C and O...
        }

Captures
--------

    captures(mol, type = Unique)
    captures_unique(mol)
    captures_all(mol)

    captures_atom(mol, atom, type = Unique)
    captures_atom_unique(mol, atom)
    captures_atom_all(mol, atom)

    captures_bond(mol, bond, type = Unique)
    captures_bond_unique(mol)
    captures_bond_all(mol)

    captures(mol, atom/bond, type = Unique)
    captures_unique(mol, atom/bond)
    captures_all(mol, atom/bond)

    Return type: std::vector<std::tuple<Atom...>>

    Example
    ~~~~~~~

        for (auto [C, N] : ctse::captures<"C-N">(mol)) {
            // Use atoms C and N
        }

*/

namespace ctse = Kitimar::CTSmarts;
