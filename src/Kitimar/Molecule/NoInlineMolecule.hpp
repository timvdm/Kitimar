#pragma once

#include "Molecule.hpp"

#include <span>

namespace Kitimar::Molecule {

    struct NoInlineAtom {};
    struct NoInlineBond {};
    struct NoInlineMolecule {};

    bool operator!=(const NoInlineAtom &a, const NoInlineAtom &b);

    // AtomList
    int num_atoms(const NoInlineMolecule&);
    NoInlineAtom get_atom(const NoInlineMolecule&, int);
    std::span<NoInlineAtom> get_atoms(const NoInlineMolecule&);
    int get_index(const NoInlineMolecule&, const NoInlineAtom&);
    NoInlineAtom null_atom(const NoInlineMolecule&);

    // BondList
    int num_bonds(const NoInlineMolecule&);
    NoInlineBond get_bond(const NoInlineMolecule&, int);
    auto get_bonds(const NoInlineMolecule&) -> std::span<NoInlineBond>;
    int get_index(const NoInlineMolecule&, const NoInlineBond&);
    NoInlineAtom get_source(const NoInlineMolecule&, const NoInlineBond&);
    NoInlineAtom get_target(const NoInlineMolecule&, const NoInlineBond&);
    NoInlineBond null_bond(const NoInlineMolecule&);

    // IncidentBondList
    int get_degree(const NoInlineMolecule&, const NoInlineAtom&);
    std::span<NoInlineBond> get_bonds(const NoInlineMolecule&, const NoInlineAtom&);

    // AdjacentAtomList
    std::span<NoInlineAtom> get_nbrs(const NoInlineMolecule&, const NoInlineAtom&);

    // ElementLayer
    int get_element(const NoInlineMolecule&, const NoInlineAtom&);

    // IsotopeLayer
    int get_isotope(const NoInlineMolecule&, const NoInlineAtom&);

    // ChargeLayer
    int get_charge(const NoInlineMolecule&, const NoInlineAtom&);

    // ValenceLayer
    int get_valence(const NoInlineMolecule&, const NoInlineAtom&);

    // BondOrderLayer
    int get_order(const NoInlineMolecule&, const NoInlineBond&);

    // ImplicitHydrogensLayer
    int get_implicit_hydrogens(const NoInlineMolecule&, const NoInlineAtom&);
    int get_total_hydrogens(const NoInlineMolecule&, const NoInlineAtom&);

    // MobileHydrogensLayer
    int get_mobile_hydrogens(const NoInlineMolecule&, const NoInlineAtom&);

    // RingLayer
    bool is_ring_atom(const NoInlineMolecule&, const NoInlineAtom&);
    bool is_ring_bond(const NoInlineMolecule&, const NoInlineBond&);

    // RingSetLayer
    bool is_in_ring_size(const NoInlineMolecule&, const NoInlineAtom&, int);
    int get_ring_count(const NoInlineMolecule&, const NoInlineAtom&, int);
    int get_ring_degree(const NoInlineMolecule&, const NoInlineAtom&, int);

    // AromaticLayer
    bool is_aromatic_atom(const NoInlineMolecule&, const NoInlineAtom&);
    bool is_aromatic_bond(const NoInlineMolecule&, const NoInlineBond&);


    static_assert(Kitimar::Molecule::Molecule<NoInlineMolecule>);

} // namespace Kitimar::Molecule
