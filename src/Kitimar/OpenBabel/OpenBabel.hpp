#pragma once

#include <Kitimar/Molecule/Molecule.hpp>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>

#include <functional>
#include <ranges>
#include <iostream>

namespace Kitimar {

    inline OpenBabel::OBMol readSmilesOpenBabel(std::string_view smiles)
    {
        OpenBabel::OBMol mol;
        OpenBabel::OBConversion conv;
        conv.SetInFormat("smi");
        conv.ReadString(&mol, smiles.data());
        return mol;
    }

    inline std::string writeSmiles(OpenBabel::OBMol &mol)
    {
        OpenBabel::OBConversion conv;
        conv.SetOutFormat("smi");
        return conv.WriteString(&mol, true);
    }

    class OpenBabelSmilesMolSource
    {
        public:
            OpenBabelSmilesMolSource(const std::string_view &filename)
            {
                m_atEnd = !m_conv.ReadFile(&m_mol, filename.data());
            }

            operator bool()
            {
                return !m_atEnd;
            }

            auto index() const noexcept
            {
                return m_index;
            }

            OpenBabel::OBConversion& obConversion()
            {
                return m_conv;
            }

            OpenBabel::OBMol read()
            {
                if (m_atEnd)
                    return {};
                auto mol = m_mol; // FIXME use shared_ptr
                m_atEnd = !m_conv.Read(&m_mol);
                ++m_index;
                return mol;
            }

            auto molecules()
            {
                return std::views::iota(0) |
                        std::views::take_while([this] (auto) {
                            if (m_atEnd)
                                m_numMolecules = m_index;
                            return !m_atEnd;
                        }) |
                        std::views::transform([this] (auto) {                            
                            return read();
                        });
            }

            void reset()
            {
                m_index = 0;
                m_atEnd = !m_conv.ReadFile(&m_mol, m_conv.GetInFilename());
            }

            auto numMolecules()
            {
                if (m_numMolecules == -1) {
                    auto index = m_index;
                    auto pos = m_conv.GetInStream()->tellg();
                    reset();

                    for (auto mol : molecules()) {
                        //m_numMolecules = std::ranges::distance(molecules());
                        //if (m_index % 1000 == 0) std::cout << m_index << std::endl;
                    }
                    m_numMolecules = m_index;

                    reset();
                    m_index = index;
                    m_conv.GetInStream()->seekg(pos);
                }
                return m_numMolecules;
            }

        private:
            OpenBabel::OBMol m_mol;
            OpenBabel::OBConversion m_conv;
            std::size_t m_index = 0;
            std::size_t m_numMolecules = -1;
            bool m_atEnd = false;
    };

} // namespace Kitimar

namespace OpenBabel {

    // AtomList

    inline auto num_atoms(const OpenBabel::OBMol &mol)
    {
        return mol.NumAtoms();
    }

    inline auto get_atoms(OpenBabel::OBMol &mol)
    {
        return std::ranges::subrange(mol.BeginAtoms(), mol.EndAtoms());
    }

    inline auto get_atom(const OpenBabel::OBMol &mol, auto index)
    {
        assert(index >= 0 && index < num_atoms(mol));
        return mol.GetAtom(index + 1);
    }

    inline auto get_index(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetIndex();
    }

    inline auto get_degree(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->GetExplicitDegree();
    }

    inline auto null_atom(const OpenBabel::OBMol &mol)
    {
        return nullptr;
    }

    // BondList

    inline auto num_bonds(const OpenBabel::OBMol &mol)
    {
        return mol.NumBonds();
    }

    inline auto get_bonds(OpenBabel::OBMol &mol)
    {
        return std::ranges::subrange(mol.BeginBonds(), mol.EndBonds());
    }

    inline auto get_bond(const OpenBabel::OBMol &mol, auto index)
    {
        assert(index >= 0 && index < num_bonds(mol));
        return mol.GetBond(index);
    }

    inline auto get_index(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetIdx();
    }

    inline auto get_source(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetBeginAtom();
    }

    inline auto get_target(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetEndAtom();
    }

    inline auto null_bond(const OpenBabel::OBMol &mol)
    {
        return nullptr;
    }

    // IncidentBondList

    inline auto get_bonds(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return std::ranges::subrange(atom->BeginBonds(), atom->EndBonds());
    }

    // AdjacentAtomList

    inline auto get_nbrs(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) { return Kitimar::Molecule::get_nbr(mol, bond, atom); });
    }

    // ElementLayer

    inline auto get_element(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetAtomicNum();
    }

    // IsotopeLayer

    inline auto get_isotope(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetIsotope();
    }

    // ChargeLayer

    inline auto get_charge(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetFormalCharge();
    }

    // ValenceLayer

    inline auto get_valence(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->GetTotalValence();
    }

    // BonderOrderLayer

    inline auto get_order(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetBondOrder();
    }

    // ImplicitHydrogenLayer

    inline auto get_implicit_hydrogens(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->GetImplicitHCount();
    }

    inline auto get_total_hydrogens(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->ExplicitHydrogenCount() + atom->GetImplicitHCount();
    }

    // RingLayer

    inline auto is_ring_atom(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->IsInRing();
    }

    inline auto is_ring_bond(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        return bond->IsInRing();
    }

    // AromaticLayer

    inline auto is_aromatic_atom(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->IsAromatic();
    }

    inline auto is_aromatic_bond(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        return bond->IsAromatic();
    }

    // RingSetLayer

    // is the atom part of a ring with the specified size
    inline auto is_in_ring_size(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom, int ringSize)
    {
        return atom->IsInRingSize(ringSize);
    }

    // number of rings the atom is part of (depends on used ring set!)
    inline auto get_ring_count(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->MemberOfRingCount();
    }

    // number of ring bonds around the atom
    inline auto get_ring_degree(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->CountRingBonds();
    }

} // namespace OpenBabel
