#pragma once

#include "Molecule.hpp"

#include <cstdint>
#include <vector>
#include <cassert>
#include <limits>

namespace Kitimar::Molecule {



    struct MockSmallAtom
    {
        uint8_t element;
        uint8_t isotope;
        int8_t charge;
        uint8_t degree;
        uint8_t implicitHydrogens;
        uint8_t explicitHydrogens;
        bool cyclic;
        bool aromatic;
    };

    struct MockSmallBond
    {
        uint32_t source;
        uint32_t target;
        uint8_t order;
        bool cyclic;
        bool aromatic;
    };

    // MockIndexMolecule

    struct MockIndexMolecule
    {
        using Atom = uint32_t;
        using Bond = uint32_t;

        std::vector<MockSmallAtom> atoms;
        std::vector<MockSmallBond> bonds;
    };

    // MockProxyMolecule

    struct MockProxyMolecule;

    struct MockProxyAtom
    {
        const MockProxyMolecule *molecule = nullptr;
        uint32_t index = std::numeric_limits<uint32_t>::max();

        bool operator==(const MockProxyAtom &other) const noexcept = default;
    };

    struct MockProxyBond
    {
        const MockProxyMolecule *molecule = nullptr;
        uint32_t index = std::numeric_limits<uint32_t>::max();
    };

    struct MockProxyMolecule
    {
        using Atom = MockProxyAtom;
        using Bond = MockProxyBond;

        std::vector<MockSmallAtom> atoms;
        std::vector<MockSmallBond> bonds;
    };

    struct MockReferenceMolecule  {};
    struct MockPointerMolecule  {};


    template<typename Mol>
    concept MockMolecule = std::same_as<Mol, MockIndexMolecule> ||
                             std::same_as<Mol, MockProxyMolecule> ||
                             std::same_as<Mol, MockReferenceMolecule> ||
                             std::same_as<Mol, MockPointerMolecule>;


    template<MockMolecule Mol>
    inline auto addMockAtom(Mol &mol, uint32_t index,
                            uint8_t element, uint8_t isotope, int8_t charge,
                            uint8_t degree, uint8_t implicitHydrogens, uint8_t explicitHydrogens,
                            bool cyclic, bool aromatic)
    {
        if constexpr (std::is_same_v<Mol, MockIndexMolecule> || std::is_same_v<Mol, MockProxyMolecule>)
            mol.atoms.emplace_back(MockSmallAtom{element, isotope, charge, degree, implicitHydrogens, explicitHydrogens, cyclic, aromatic});
    }

    template<MockMolecule Mol>
    inline auto addMockBond(Mol &mol, uint32_t index,
                            uint32_t source, uint32_t target,uint8_t order,
                            bool cyclic, bool aromatic)
    {
        if constexpr (std::is_same_v<Mol, MockIndexMolecule> || std::is_same_v<Mol, MockProxyMolecule>)
            mol.bonds.emplace_back(MockSmallBond{source, target, order, cyclic, aromatic});
    }


    //
    // Molecule
    //

    inline uint32_t num_atoms(const MockMolecule auto &mol) noexcept
    {
        return mol.atoms.size();
    }

    inline uint32_t num_bonds(const MockMolecule auto &mol) noexcept
    {
        return mol.bonds.size();
    }

    template<MockMolecule Mol>
    inline Mol::Atom get_atom(const Mol &mol, uint32_t index) noexcept
    {
        assert(index < num_atoms(mol));
        if constexpr (std::is_same_v<Mol, MockIndexMolecule>)
            return index;
        else if constexpr (std::is_same_v<Mol, MockProxyMolecule>)
            return MockProxyAtom{&mol, index};
    }

    template<MockMolecule Mol>
    inline Mol::Bond get_bond(const Mol &mol, uint32_t index) noexcept
    {
        assert(index < num_bonds(mol));
        if constexpr (std::is_same_v<Mol, MockIndexMolecule>)
            return index;
        else if constexpr (std::is_same_v<Mol, MockProxyMolecule>)
            return MockProxyBond{&mol, index};
    }

    inline auto get_atoms(const MockMolecule auto &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_atoms(mol))>(0), num_atoms(mol))
                | std::views::transform([&mol] (auto index) { return get_atom(mol, index); });
    }


    inline auto get_bonds(const MockMolecule auto &mol) noexcept
    {
        return std::views::iota(static_cast<decltype(num_bonds(mol))>(0), num_bonds(mol))
                | std::views::transform([&mol] (auto index) { return get_bond(mol, index); });
    }

    template<MockMolecule Mol>
    inline auto get_index(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        if constexpr (std::is_same_v<Mol, MockIndexMolecule>)
            return atom;
        else if constexpr (std::is_same_v<Mol, MockProxyMolecule>)
            return atom.index;
    }

    template<MockMolecule Mol>
    inline auto get_index(const Mol &mol, typename Mol::Bond bond) noexcept
    {
        if constexpr (std::is_same_v<Mol, MockIndexMolecule>)
            return bond;
        else if constexpr (std::is_same_v<Mol, MockProxyMolecule>)
            return bond.index;
    }

    inline auto get_index(const MockIndexMolecule &mol, uint32_t index) noexcept
    {
        return index;
    }


    //
    // Atom
    //

    template<MockMolecule Mol>
    inline auto get_element(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].element;
    }

    template<MockMolecule Mol>
    inline auto get_isotope(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].isotope;
    }

    template<MockMolecule Mol>
    inline auto get_charge(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].charge;
    }

    template<MockMolecule Mol>
    inline auto get_degree(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].degree;
    }

    template<MockMolecule Mol>
    inline auto get_implicit_hydrogens(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].implicitHydrogens;
    }

    template<MockMolecule Mol>
    inline auto get_total_hydrogens(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].explicitHydrogens + mol.atoms[get_index(mol, atom)].implicitHydrogens;
    }

    template<MockMolecule Mol>
    inline auto is_ring_atom(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].cyclic;
    }

    template<MockMolecule Mol>
    inline auto is_aromatic_atom(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return mol.atoms[get_index(mol, atom)].aromatic;
    }

    template<MockMolecule Mol>
    inline Mol::Atom get_source(const Mol &mol, typename Mol::Bond bond) noexcept;
    template<MockMolecule Mol>
    inline Mol::Atom get_target(const Mol &mol, typename Mol::Bond bond) noexcept;

    template<MockMolecule Mol>
    inline auto get_bonds(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        //return get_bonds(mol) | std::views::filter([&mol, atom] (auto i) {
        //    auto bond = get_bond(mol, i);
        return get_bonds(mol) | std::views::filter([&mol, atom] (auto bond) {
            return atom == get_source(mol, bond) || atom == get_target(mol, bond);
        });
    }

    template<MockMolecule Mol>
    inline auto get_nbrs(const Mol &mol, typename Mol::Atom atom) noexcept
    {
        assert(get_index(mol, atom) < num_atoms(mol));
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) {
            return Kitimar::Molecule::get_nbr(mol, bond, atom);
        });
    }

    template<MockMolecule Mol>
    inline Mol::Atom null_atom(const Mol &mol) noexcept
    {
        if constexpr (std::is_same_v<Mol, MockIndexMolecule>)
            return std::numeric_limits<uint32_t>::max();
        if constexpr (std::is_same_v<Mol, MockProxyMolecule>)
            return {};
    }

    //
    // Bond
    //

    template<MockMolecule Mol>
    inline Mol::Atom get_source(const Mol &mol, typename Mol::Bond bond) noexcept
    {
        assert(get_index(mol, bond) < num_bonds(mol));
        return get_atom(mol, mol.bonds[get_index(mol, bond)].source);
    }

    template<MockMolecule Mol>
    inline Mol::Atom get_target(const Mol &mol, typename Mol::Bond bond) noexcept
    {
        assert(get_index(mol, bond) < num_bonds(mol));
        return get_atom(mol, mol.bonds[get_index(mol, bond)].target);
    }

    template<MockMolecule Mol>
    inline auto get_order(const Mol &mol, typename Mol::Bond bond) noexcept
    {
        assert(get_index(mol, bond) < num_bonds(mol));
        return mol.bonds[get_index(mol, bond)].order;
    }

    template<MockMolecule Mol>
    inline auto is_ring_bond(const Mol &mol, typename Mol::Bond bond) noexcept
    {
        assert(get_index(mol, bond) < num_bonds(mol));
        return mol.bonds[get_index(mol, bond)].cyclic;
    }

    template<MockMolecule Mol>
    inline auto is_aromatic_bond(const Mol &mol, typename Mol::Bond bond) noexcept
    {
        assert(get_index(mol, bond) < num_bonds(mol));
        return mol.bonds[get_index(mol, bond)].aromatic;
    }

    template<MockMolecule Mol>
    inline Mol::Bond null_bond(const Mol &mol) noexcept
    {
        if constexpr (std::is_same_v<Mol, MockIndexMolecule>)
            return std::numeric_limits<uint32_t>::max();
        if constexpr (std::is_same_v<Mol, MockProxyMolecule>)
            return {};
    }

} // namespace Kitimar::Molecule
