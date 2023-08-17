#pragma once

#include "API/Match.hpp"
#include "API/Count.hpp"
#include "API/Map.hpp"
#include "API/Maps.hpp"
#include "API/Capture.hpp"
#include "API/Captures.hpp"

namespace Kitimar::CTSmarts {

    /*
     * Match
     *
     * bool match(mol)
     * bool match_atom(mol, atom)
     * bool match_bond(mol, bond)
     * bool match(mol, atom/bond)
     *
     * Count
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
     * Map
     *
     * (bool, Map) map(mol)
     * (bool, Map) Map map_atom(mol, atom)
     * (bool, Map) Map map_bond(mol, bond)
     * (bool, Map) Map map(mol, atom/bond)
     *
     * Maps
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
     * Capture
     *
     * (bool, Atom...) capture(mol)
     * (bool, Atom...) capture_atom(mol, atom)
     * (bool, Atom...) capture_bond(mol, bond)
     * (bool, Atom...) capture(mol, atom/bond)
     *
     * Captures
     *
     * Maps captures(mol, type = Unique)
     * Maps captures_unique(mol)
     * Maps captures_all(mol)
     *
     * Maps captures_atom(mol, atom, type = Unique)
     * Maps captures_atom_unique(mol, atom)
     * Maps captures_atom_all(mol, atom)
     *
     * Maps captures_bond(mol, bond, type = Unique)
     * Maps captures_bond_unique(mol)
     * Maps captures_bond_all(mol)
     *
     * Maps captures(mol, atom/bond, type = Unique)
     * Maps captures_unique(mol, atom/bond)
     * Maps captures_all(mol, atom/bond)
     *
     */

} // namespace Kitimar::CTSmarts

namespace ctse = Kitimar::CTSmarts;
