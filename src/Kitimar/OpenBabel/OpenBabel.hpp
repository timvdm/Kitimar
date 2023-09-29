#pragma once

#include "../Molecule/Molecule.hpp"
#include "../Molecule/Toolkit.hpp"
#include "../Util/Util.hpp"

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/parsmart.h>

#include <functional>
#include <ranges>
#include <iostream>

namespace Kitimar {

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

            std::shared_ptr<OpenBabel::OBMol> read()
            {
                if (m_atEnd)
                    return {};
                auto mol = std::make_shared<OpenBabel::OBMol>(m_mol);
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
                if (m_numMolecules == static_cast<std::size_t>(-1)) {
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

    namespace Toolkit {

        static constexpr inline auto openbabel = toolkitId("OpenBabel");

        template<>
        inline auto name<openbabel>()
        {
            return "OpenBabel";
        }

        template<>
        inline auto readSmiles<openbabel>(std::string_view smiles)
        {
            auto mol = std::make_shared<OpenBabel::OBMol >();
            OpenBabel::OBConversion conv;
            conv.SetInFormat("smi");
            conv.ReadString(mol.get(), smiles.data());
            return mol;
        }

        template<>
        inline auto writeSmiles<openbabel, OpenBabel::OBMol>(const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBConversion conv;
            conv.SetOutFormat("smi");
            return conv.WriteString(const_cast<OpenBabel::OBMol*>(&mol), true);
        }

        template<>
        inline auto smilesMolSource<openbabel>(std::string_view filename)
        {
            return OpenBabelSmilesMolSource{filename};
        }

        namespace impl {

            inline auto openbabelSmartsAll(std::string_view SMARTS, const OpenBabel::OBMol &mol, bool unique)
            {
                OpenBabel::OBSmartsPattern smarts;
                smarts.Init(std::string(SMARTS));
                smarts.Match(const_cast<OpenBabel::OBMol&>(mol));
                return unique ? smarts.GetUMapList() : smarts.GetMapList();
            }

        } // namespace impl

        template<>
        inline auto match<openbabel, OpenBabel::OBMol>(std::string_view SMARTS, const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBSmartsPattern smarts;
            smarts.Init(std::string(SMARTS));
            return smarts.HasMatch(const_cast<OpenBabel::OBMol&>(mol));
        }

        template<>
        inline auto map<openbabel, OpenBabel::OBMol>(std::string_view SMARTS, const OpenBabel::OBMol &mol)
        {
            OpenBabel::OBSmartsPattern smarts;
            smarts.Init(std::string(SMARTS));
            auto found = smarts.Match(const_cast<OpenBabel::OBMol&>(mol));
            if (found) {
                auto maps = smarts.GetMapList();
                assert(maps.size() == 1);
                return std::make_tuple(found, maps[0]);
            }
            return std::make_tuple(found, std::vector<int>{});
        }

        template<>
        inline auto count_unique<openbabel, OpenBabel::OBMol>(std::string_view SMARTS, const OpenBabel::OBMol &mol)
        {
            return impl::openbabelSmartsAll(SMARTS, mol, true).size();
        }

        template<>
        inline auto count_all<openbabel, OpenBabel::OBMol>(std::string_view SMARTS, const OpenBabel::OBMol &mol)
        {
            return impl::openbabelSmartsAll(SMARTS, mol, false).size();
        }

        template<>
        inline auto maps_unique<openbabel, OpenBabel::OBMol>(std::string_view SMARTS, const OpenBabel::OBMol &mol)
        {
            return impl::openbabelSmartsAll(SMARTS, mol, true);
        }

        template<>
        inline auto maps_all<openbabel, OpenBabel::OBMol>(std::string_view SMARTS, const OpenBabel::OBMol &mol)
        {
            return impl::openbabelSmartsAll(SMARTS, mol, false);
        }


    } // namespace Toolkit

} // namespace Kitimar

namespace OpenBabel {

    // AtomList

    inline auto num_atoms(const OpenBabel::OBMol &mol)
    {
        return mol.NumAtoms();
    }

    inline auto get_atoms(const OpenBabel::OBMol &mol)
    {
        auto &m = const_cast<OpenBabel::OBMol&>(mol);
        return std::ranges::subrange(m.BeginAtoms(), m.EndAtoms());
    }

    inline auto get_atom(const OpenBabel::OBMol &mol, auto index)
    {
        assert(index != static_cast<decltype(index)>(-1));
        assert(static_cast<unsigned int>(index) < num_atoms(mol));
        return mol.GetAtom(index + 1);
    }

    inline auto get_index(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetIndex();
    }

    inline auto get_degree(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->GetExplicitDegree();
    }

    inline OpenBabel::OBAtom* null_atom(const OpenBabel::OBMol&)
    {
        return nullptr;
    }

    // BondList

    inline auto num_bonds(const OpenBabel::OBMol &mol)
    {
        return mol.NumBonds();
    }

    inline auto get_bonds(const OpenBabel::OBMol &mol)
    {
        auto &m = const_cast<OpenBabel::OBMol&>(mol);
        return std::ranges::subrange(m.BeginBonds(), m.EndBonds());
    }

    inline auto get_bond(const OpenBabel::OBMol &mol, auto index)
    {
        assert(index != static_cast<decltype(index)>(-1));
        assert(static_cast<unsigned int>(index) < num_bonds(mol));
        return mol.GetBond(index);
    }

    inline auto get_index(const OpenBabel::OBMol&, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetIdx();
    }

    inline auto get_source(const OpenBabel::OBMol&, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetBeginAtom();
    }

    inline auto get_target(const OpenBabel::OBMol&, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetEndAtom();
    }

    inline OpenBabel::OBBond* null_bond(const OpenBabel::OBMol&)
    {
        return nullptr;
    }

    // IncidentBondList

    inline auto get_bonds(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
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

    inline auto get_element(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetAtomicNum();
    }

    // IsotopeLayer

    inline auto get_isotope(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetIsotope();
    }

    // ChargeLayer

    inline auto get_charge(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        assert(atom);
        return atom->GetFormalCharge();
    }

    // ValenceLayer

    inline auto get_valence(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->GetTotalValence();
    }

    // BonderOrderLayer

    inline auto get_order(const OpenBabel::OBMol&, OpenBabel::OBBond *bond)
    {
        assert(bond);
        return bond->GetBondOrder();
    }

    // ImplicitHydrogenLayer

    inline auto get_implicit_hydrogens(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->GetImplicitHCount();
    }

    inline auto get_total_hydrogens(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        // Avoid function call overhead to slow atom->ExplicitHydrogenCount()
        auto count = 0;
        for (auto nbr : get_nbrs(mol, atom))
            if (get_element(mol, nbr) == 1)
                ++count;
        return count + atom->GetImplicitHCount();
    }

    // RingLayer

    inline auto is_ring_atom(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->IsInRing();
    }

    inline auto is_ring_bond(const OpenBabel::OBMol&, OpenBabel::OBBond *bond)
    {
        return bond->IsInRing();
    }

    // AromaticLayer

    inline auto is_aromatic_atom(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->IsAromatic();
    }

    inline auto is_aromatic_bond(const OpenBabel::OBMol&, OpenBabel::OBBond *bond)
    {
        return bond->IsAromatic();
    }

    // RingSetLayer

    // is the atom part of a ring with the specified size
    inline auto is_in_ring_size(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom, int ringSize)
    {
        return atom->IsInRingSize(ringSize);
    }

    // number of rings the atom is part of (depends on used ring set!)
    inline auto get_ring_count(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->MemberOfRingCount();
    }

    // number of ring bonds around the atom
    inline auto get_ring_degree(const OpenBabel::OBMol&, OpenBabel::OBAtom *atom)
    {
        return atom->CountRingBonds();
    }

} // namespace OpenBabel
