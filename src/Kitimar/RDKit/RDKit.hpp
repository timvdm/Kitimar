#pragma once

#include <Kitimar/Molecule/Molecule.hpp>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>

#include <ranges>

namespace Kitimar {

    inline std::shared_ptr<RDKit::ROMol> readSmilesRDKit(std::string_view smiles)
    {
        return std::shared_ptr<RDKit::ROMol>{RDKit::SmilesToMol(smiles.data())};
    }

    class RDKitSmilesMolSource
    {
        public:
            RDKitSmilesMolSource(const std::string_view &filename) : m_supplier{filename.data(), "\t", 0, 1, false}
            {
            }

            operator bool()
            {
                return m_atEnd;
            }

            auto index() const noexcept
            {
                return m_index;
            }

            RDKit::SmilesMolSupplier& supplier()
            {
                return m_supplier;
            }

            std::shared_ptr<RDKit::ROMol> read()
            {
                if (m_supplier.atEnd()) {
                    m_atEnd = true;
                    return nullptr;
                }
                ++m_index;

                auto mol = std::shared_ptr<RDKit::ROMol>{m_supplier.next()};
                m_atEnd = m_supplier.atEnd();
                return mol;

            }

            auto molecules()
            {
                return std::views::iota(0) |
                        std::views::take_while([this] (auto i) {
                            return !m_atEnd;
                        }) |
                        std::views::transform([this] (auto i) {
                            return read();
                        });
            }

        private:
            RDKit::SmilesMolSupplier m_supplier;
            std::size_t m_index = 0;
            bool m_atEnd = false;
    };

    class RDKitPickleMolSource
    {
        public:
            RDKitPickleMolSource(const std::string_view &filename) : m_ifs{filename.data(), std::ios_base::binary}
            {                
                m_atEnd = m_ifs.eof();
                if (!m_atEnd) {
                    m_mol.reset(new RDKit::ROMol);
                    RDKit::MolPickler::molFromPickle(m_ifs, *m_mol);
                }
            }

            operator bool()
            {
                return m_atEnd;
            }

            auto index() const noexcept
            {
                return m_index;
            }

            auto& ifs()
            {
                return m_ifs;
            }

            std::shared_ptr<RDKit::ROMol> read()
            {
                if (m_atEnd)
                    return nullptr;

                auto mol = m_mol;
                try {
                    m_mol.reset(new RDKit::ROMol);
                    RDKit::MolPickler::molFromPickle(m_ifs, *m_mol);
                } catch( RDKit::MolPicklerException &e ) {
                    m_atEnd = true;
                }

                ++m_index;
                return mol;

            }

            auto molecules()
            {
                return std::views::iota(0) |
                        std::views::take_while([this] (auto i) {
                            return !m_atEnd;
                        }) |
                        std::views::transform([this] (auto i) {
                            return read();
                        });
            }

        private:
            std::shared_ptr<RDKit::ROMol> m_mol;
            std::ifstream m_ifs;
            std::size_t m_index = 0;
            bool m_atEnd = false;
    };


} // namespace Kitimar


namespace RDKit {

    // AtomList

    inline auto num_atoms(const RDKit::ROMol &mol)
    {
        return mol.getNumAtoms();
    }

    inline auto get_atoms(const RDKit::ROMol &mol)
    {
        return std::views::iota(static_cast<decltype(mol.getNumAtoms())>(0), mol.getNumAtoms()) |
                std::views::transform([&mol] (auto idx) {
                    return mol.getAtomWithIdx(idx);
                });
    }

    inline auto get_atom(const RDKit::ROMol &mol, auto index)
    {
        assert(index >= 0 && index < num_atoms(mol));
        return mol.getAtomWithIdx(index);
    }

    inline auto get_index(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getIdx();
    }

    inline auto get_degree(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        return atom->getDegree();
    }

    inline RDKit::Atom* null_atom(const RDKit::ROMol &mol)
    {
        return nullptr;
    }

    // BondList

    inline auto num_bonds(const RDKit::ROMol &mol)
    {
        return mol.getNumBonds();
    }

    inline auto get_bonds(const RDKit::ROMol &mol)
    {
        return std::views::iota(static_cast<decltype(mol.getNumBonds())>(0), mol.getNumBonds()) |
                std::views::transform([&mol] (auto idx) {
                    return mol.getBondWithIdx(idx);
                });
    }

    inline auto get_bond(const RDKit::ROMol &mol, auto index)
    {
        assert(index >= 0 && index < num_bonds(mol));
        return mol.getBondWithIdx(index);
    }

    inline auto get_index(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getIdx();
    }

    inline auto get_source(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getBeginAtom();
    }

    inline auto get_target(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getEndAtom();
    }

    inline RDKit::Bond* null_bond(const RDKit::ROMol &mol)
    {
        return nullptr;
    }


    // IncidentBondList

    inline auto get_bonds(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        return boost::make_iterator_range(mol.getAtomBonds(atom)) |
                std::views::transform([&mol] (auto idx) {
                    return mol[idx];
                });
    }

    // AdjacentAtomList


    inline auto get_nbrs(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return get_bonds(mol, atom) |
                std::views::transform([atom] (const RDKit::Bond *bond) { return bond->getOtherAtom(atom); });
    }

    // ElementLayer

    inline auto get_element(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getAtomicNum();
    }

    // IsotopeLayer

    inline auto get_isotope(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getIsotope();
    }

    // ChargeLayer

    inline auto get_charge(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getFormalCharge();
    }

    // BonderOrderLayer

    inline auto get_order(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        assert(bond);
        switch (bond->getBondType()) {
            case RDKit::Bond::SINGLE: return 1;
            case RDKit::Bond::DOUBLE: return 2;
            case RDKit::Bond::TRIPLE: return 3;
            case RDKit::Bond::QUADRUPLE: return 4;
            default:
                return 1;
        }
    }

    // ImplicitHydrogenLayer

    inline auto get_implicit_hydrogens(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        return atom->getImplicitValence();
    }

    // AromaticLayer

    inline auto is_aromatic_atom(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        return atom->getIsAromatic();
    }

    inline auto is_aromatic_bond(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        return bond->getIsAromatic();
    }

    inline auto is_ring_bond(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        if (!mol.getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(mol);
        return mol.getRingInfo()->numBondRings(bond->getIdx()) > 0;
    }

}
