#pragma once

#include "Molecule.hpp"

#include <cstdint>
#include <vector>
#include <cassert>

namespace Kitimar::Molecule {

    struct MockAtom
    {
        uint8_t element;
        uint8_t isotope;
        int8_t charge;
        uint8_t degree;
        uint8_t implicitHydrogens;
        bool cyclic;
        bool aromatic;
    };

    struct MockBond
    {
        uint32_t source;
        uint32_t target;
        uint8_t order;
        bool cyclic;
        bool aromatic;
    };


    struct MockMolecule
    {
        std::vector<MockAtom> atoms;
        std::vector<MockBond> bonds;
    };

    //
    // Molecule
    //

    inline auto num_atoms(const MockMolecule &mol) noexcept
    {
        return mol.atoms.size();
    }

    inline auto num_bonds(const MockMolecule &mol) noexcept
    {
        return mol.atoms.size();
    }

    inline auto get_atoms(const MockMolecule &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_atoms(mol))>(0), num_atoms(mol));
    }

    inline auto get_bonds(const MockMolecule &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_bonds(mol))>(0), num_bonds(mol));
    }

    inline auto get_atom(const MockMolecule &mol, uint32_t index) noexcept
    {
        assert(index < num_atoms(mol));
        return index;
    }

    inline auto get_bond(const MockMolecule &mol, uint32_t index) noexcept
    {
        assert(index < num_bonds(mol));
        return index;
    }

    inline auto get_index(const MockMolecule &mol, uint32_t index) noexcept
    {
        return index;
    }

    //
    // Atom
    //

    inline auto get_element(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].element;
    }

    inline auto get_isotope(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].isotope;
    }

    inline auto get_charge(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].charge;
    }

    inline auto get_degree(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].degree;
    }

    inline auto get_implicit_hydrogens(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].implicitHydrogens;
    }

    inline auto is_cyclic_atom(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].cyclic;
    }

    inline auto is_aromatic_atom(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return mol.atoms[atom].aromatic;
    }

    inline uint32_t get_source(const MockMolecule &mol, uint32_t bond) noexcept;
    inline uint32_t get_target(const MockMolecule &mol, uint32_t bond) noexcept;

    inline auto get_bonds(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return get_bonds(mol) | std::views::filter([&mol, atom] (auto i) {
            auto bond = get_bond(mol, i);
            return atom == get_source(mol, bond) || atom == get_target(mol, bond);
        });
    }

    inline auto get_nbrs(const MockMolecule &mol, uint32_t atom) noexcept
    {
        assert(atom < num_atoms(mol));
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
            return Kitimar::Molecule::get_nbr(mol, bond, atom);
        });
    }

    inline auto null_atom(const MockMolecule &mol) noexcept
    {
        return -1;
    }

    //
    // Bond
    //

    inline uint32_t get_source(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].source;
    }

    inline uint32_t get_target(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].target;
    }

    inline auto get_order(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].order;
    }

    inline auto is_cyclic_bond(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].cyclic;
    }

    inline auto is_aromatic_bond(const MockMolecule &mol, uint32_t bond) noexcept
    {
        assert(bond < num_bonds(mol));
        return mol.bonds[bond].aromatic;
    }

    inline auto null_bond(const MockMolecule &mol) noexcept
    {
        return -1;
    }

} // namespace Kitimar::Molecule
