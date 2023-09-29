#pragma once

#include "../Molecule/Molecule.hpp"
#include "../Molecule/Toolkit.hpp"

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <ranges>

namespace Kitimar {

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
                        std::views::take_while([this] (auto) {
                            return !m_atEnd;
                        }) |
                        std::views::transform([this] (auto) {
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
                        std::views::take_while([this] (auto) {
                            return !m_atEnd;
                        }) |
                        std::views::transform([this] (auto) {
                            return read();
                        });
            }

        private:
            std::shared_ptr<RDKit::ROMol> m_mol;
            std::ifstream m_ifs;
            std::size_t m_index = 0;
            bool m_atEnd = false;
    };

    namespace Toolkit {

        static constexpr inline auto rdkit = toolkitId("RDKit");

        template<>
        inline auto name<rdkit>()
        {
            return "RDKit";
        }

        template<>
        inline auto readSmiles<rdkit>(std::string_view smiles)
        {
            return std::shared_ptr<RDKit::ROMol>{RDKit::SmilesToMol(smiles.data())};
        }

        template<>
        inline auto writeSmiles<rdkit, RDKit::ROMol>(const RDKit::ROMol &mol)
        {
            return RDKit::MolToSmiles(mol);
        }

        template<>
        inline auto smilesMolSource<rdkit>(std::string_view filename)
        {
            return RDKitSmilesMolSource{filename};
        }

        namespace impl {

            inline auto rdkitSmartsSingle(std::string_view SMARTS, const RDKit::ROMol &mol)
            {
                std::unique_ptr<RDKit::RWMol> smarts{RDKit::SmartsToMol(std::string(SMARTS))};
                assert(smarts);
                RDKit::MatchVectType res;
                auto found = RDKit::SubstructMatch(mol, *smarts, res);
                return std::make_tuple(found, res);
            }

            inline auto rdkitSmartsAll(std::string_view SMARTS, const RDKit::ROMol &mol, bool unique)
            {
                std::unique_ptr<RDKit::RWMol> smarts{RDKit::SmartsToMol(std::string(SMARTS))};
                assert(smarts);
                std::vector<RDKit::MatchVectType> res;
                RDKit::SubstructMatch(mol, *smarts, res, unique);
                return res;
            }

        } // namespace impl

        template<>
        inline auto match<rdkit, RDKit::ROMol>(std::string_view SMARTS, const RDKit::ROMol &mol)
        {
            return std::get<0>(impl::rdkitSmartsSingle(SMARTS, mol));
        }

        template<>
        inline auto map<rdkit, RDKit::ROMol>(std::string_view SMARTS, const RDKit::ROMol &mol)
        {
            return impl::rdkitSmartsSingle(SMARTS, mol);
        }

        template<>
        inline auto count_unique<rdkit, RDKit::ROMol>(std::string_view SMARTS, const RDKit::ROMol &mol)
        {
            return impl::rdkitSmartsAll(SMARTS, mol, true).size();
        }

        template<>
        inline auto count_all<rdkit, RDKit::ROMol>(std::string_view SMARTS, const RDKit::ROMol &mol)
        {
            return impl::rdkitSmartsAll(SMARTS, mol, false).size();
        }

        template<>
        inline auto maps_unique<rdkit, RDKit::ROMol>(std::string_view SMARTS, const RDKit::ROMol &mol)
        {
            return impl::rdkitSmartsAll(SMARTS, mol, true);
        }

        template<>
        inline auto maps_all<rdkit, RDKit::ROMol>(std::string_view SMARTS, const RDKit::ROMol &mol)
        {
            return impl::rdkitSmartsAll(SMARTS, mol, false);
        }



    } // namespace Toolkit


} // namespace Kitimar


namespace RDKit {

    // AtomList

    inline auto num_atoms(const RDKit::ROMol &mol)
    {
        return mol.getNumAtoms();
    }

    inline auto get_atoms(const RDKit::ROMol &mol)
    {
        const auto &g = mol.getTopology();
        auto [begin, end] = boost::vertices(g);
        return std::ranges::subrange(begin, end) | std::views::transform([&g] (const auto &v) -> const RDKit::Atom* {
            return g[v];
        });
    }

    inline auto get_atom(const RDKit::ROMol &mol, auto index)
    {
        assert(index != static_cast<decltype(index)>(-1));
        assert(static_cast<unsigned int>(index) < num_atoms(mol));
        return mol.getAtomWithIdx(index);
    }

    inline auto get_index(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getIdx();
    }

    inline auto get_degree(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        return atom->getDegree();
    }

    inline const RDKit::Atom* null_atom(const RDKit::ROMol&)
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
        const auto &g = mol.getTopology();
        auto [begin, end] = boost::edges(g);
        return std::ranges::subrange(begin, end) | std::views::transform([&g] (const auto &e) -> const RDKit::Bond* {
            return g[e];
        });
    }

    inline auto get_bond(const RDKit::ROMol &mol, auto index)
    {
        assert(index != static_cast<decltype(index)>(-1));
        assert(static_cast<unsigned int>(index) < num_bonds(mol));
        return mol.getBondWithIdx(index);
    }

    inline auto get_index(const RDKit::ROMol&, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getIdx();
    }

    inline auto get_source(const RDKit::ROMol&, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getBeginAtom();
    }

    inline auto get_target(const RDKit::ROMol&, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getEndAtom();
    }

    inline const RDKit::Bond* null_bond(const RDKit::ROMol&)
    {
        return nullptr;
    }


    // IncidentBondList

    inline auto get_bonds(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        auto [begin, end] = mol.getAtomBonds(atom);
        return std::ranges::subrange(begin, end) | std::views::transform([&mol] (const auto &e) -> const RDKit::Bond* {
            return mol[e];
        });
    }

    inline auto get_nbr(const RDKit::ROMol&, const RDKit::Bond *bond, const RDKit::Atom *atom) noexcept
    {
        return atom->getIdx() == bond->getBeginAtomIdx() ? bond->getEndAtom() : bond->getBeginAtom();
    }

    // AdjacentAtomList

    inline auto get_nbrs(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (const RDKit::Bond *bond) -> const RDKit::Atom* {
            return get_nbr(mol, bond, atom);
        });
    }

    // ElementLayer

    inline auto get_element(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getAtomicNum();
    }

    // IsotopeLayer

    inline auto get_isotope(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getIsotope();
    }

    // ChargeLayer

    inline auto get_charge(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getFormalCharge();
    }

    // ValenceLayer

    inline auto get_valence(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        return atom->getTotalValence();
    }

    // BonderOrderLayer

    inline auto get_order(const RDKit::ROMol&, const RDKit::Bond *bond)
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

    inline auto get_implicit_hydrogens(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        return atom->getNumExplicitHs() + atom->getImplicitValence();
    }

    inline auto get_total_hydrogens(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        auto count = 0;
        for (auto nbr : get_nbrs(mol, atom))
            if (get_element(mol, nbr) == 1)
                ++count;
        return count + get_implicit_hydrogens(mol, atom);
    }

    // RingLayer

    inline auto is_ring_atom(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        if (!mol.getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(mol);
        return mol.getRingInfo()->numAtomRings(atom->getIdx()) > 0;
    }

    inline auto is_ring_bond(const RDKit::ROMol &mol, const RDKit::Bond *bond)
    {
        if (!mol.getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(mol);
        return mol.getRingInfo()->numBondRings(bond->getIdx()) > 0;
    }

    // AromaticLayer

    inline auto is_aromatic_atom(const RDKit::ROMol&, const RDKit::Atom *atom)
    {
        return atom->getIsAromatic();
    }

    inline auto is_aromatic_bond(const RDKit::ROMol&, const RDKit::Bond *bond)
    {
        switch (bond->getBondType()) {
            case RDKit::Bond::SINGLE:
            case RDKit::Bond::DOUBLE:
            case RDKit::Bond::AROMATIC:
                return bond->getIsAromatic();
            default:
                return false;
        }
    }

    // RingSetLayer

    // is the atom part of a ring with the specified size
    inline auto is_in_ring_size(const RDKit::ROMol &mol, const RDKit::Atom *atom, unsigned int ringSize)
    {
        if (!mol.getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(mol);
        return mol.getRingInfo()->minAtomRingSize(atom->getIdx()) == ringSize;
    }

    // number of rings the atom is part of (depends on used ring set!)
    inline auto get_ring_count(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        if (!mol.getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(mol);
        return mol.getRingInfo()->numAtomRings(atom->getIdx());
    }

    // number of ring bonds around the atom
    inline auto get_ring_degree(const RDKit::ROMol &mol, const RDKit::Atom *atom)
    {
        if (!mol.getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(mol);
        return queryAtomRingBondCount(atom);
    }

}
