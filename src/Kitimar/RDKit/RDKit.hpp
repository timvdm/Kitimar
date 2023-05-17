#pragma once

#include <Kitimar/Molecule/Molecule.hpp>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/MolOps.h>

#include <ranges>

namespace Kitimar {

    inline void readMoleculesRDKit(const std::string_view &filename,
                                   std::function<bool (std::shared_ptr<RDKit::ROMol>)> callback)
    {

        RDKit::SmilesMolSupplier supplier{filename.data()};
        while (!supplier.atEnd()) {
            std::shared_ptr<RDKit::ROMol> mol{supplier.next()};
            if (!callback(mol))
                break;
        }
    }

    inline std::shared_ptr<RDKit::ROMol> readSmilesRDKit(std::string_view smiles)
    {
        return std::shared_ptr<RDKit::ROMol>{RDKit::SmilesToMol(smiles.data())};
    }

} // namespace Kitimar


namespace RDKit {

    // AtomList

    inline auto num_atoms(const RDKit::ROMol *mol)
    {
        return mol->getNumAtoms();
    }

    inline auto get_atoms(const RDKit::ROMol *mol)
    {
        return std::views::iota(static_cast<decltype(mol->getNumAtoms())>(0), mol->getNumAtoms()) |
                std::views::transform([mol] (auto idx) {
                    return mol->getAtomWithIdx(idx);
                });
    }

    inline auto get_atom(const RDKit::ROMol *mol, auto index)
    {
        assert(index >= 0 && index < num_atoms(mol));
        return mol->getAtomWithIdx(index);
    }

    inline auto get_index(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getIdx();
    }

    inline auto get_degree(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        return atom->getDegree();
    }

    // BondList

    inline auto num_bonds(const RDKit::ROMol *mol)
    {
        return mol->getNumBonds();
    }

    inline auto get_bonds(const RDKit::ROMol *mol)
    {
        return std::views::iota(static_cast<decltype(mol->getNumBonds())>(0), mol->getNumBonds()) |
                std::views::transform([mol] (auto idx) {
                    return mol->getBondWithIdx(idx);
                });
    }

    inline auto get_bond(const RDKit::ROMol *mol, auto index)
    {
        assert(index >= 0 && index < num_bonds(mol));
        return mol->getBondWithIdx(index);
    }

    inline auto get_index(const RDKit::ROMol *mol, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getIdx();
    }

    inline auto get_source(const RDKit::ROMol *mol, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getBeginAtom();
    }

    inline auto get_target(const RDKit::ROMol *mol, const RDKit::Bond *bond)
    {
        assert(bond);
        return bond->getEndAtom();
    }












    // IncidentBondList

    inline auto get_bonds(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        return boost::make_iterator_range(mol->getAtomBonds(atom)) |
                std::views::transform([mol] (auto idx) {
                    return (*mol)[idx];
                });
    }

    // AdjacentAtomList


    inline auto get_nbrs(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return get_bonds(mol, atom) |
                std::views::transform([atom] (const RDKit::Bond *bond) { return bond->getOtherAtom(atom); });
    }

    // ElementLayer

    inline auto get_element(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getAtomicNum();
    }

    // IsotopeLayer

    inline auto get_isotope(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getIsotope();
    }

    // ChargeLayer

    inline auto get_charge(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        assert(atom);
        return atom->getFormalCharge();
    }

    // BonderOrderLayer

    inline auto get_order(const RDKit::ROMol *mol, const RDKit::Bond *bond)
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

    inline auto get_implicit_hydrogens(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        return atom->getImplicitValence();
    }

    // AromaticLayer

    inline auto is_aromatic_atom(const RDKit::ROMol *mol, const RDKit::Atom *atom)
    {
        return atom->getIsAromatic();
    }

    inline auto is_aromatic_bond(const RDKit::ROMol *mol, const RDKit::Bond *bond)
    {
        return bond->getIsAromatic();
    }

    inline auto is_cyclic_bond(const RDKit::ROMol *mol, const RDKit::Bond *bond)
    {
        if (!mol->getRingInfo()->isInitialized())
            RDKit::MolOps::findSSSR(*mol);
        return mol->getRingInfo()->numBondRings(bond->getIdx()) > 0;
    }

}
