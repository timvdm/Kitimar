#pragma once

#include <Kitimar/Molecule/Molecule.hpp>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>

#include <functional>
#include <ranges>

namespace Kitimar {

    namespace detail {

        inline OpenBabel::OBConversion makeOBConversion(const std::string &filename)
        {
            OpenBabel::OBConversion conv;
            auto *format = conv.FormatFromExt(filename);
            if (!format)
                throw std::runtime_error("Could not determine format");
            if (!conv.SetInFormat(format))
                throw std::runtime_error("Could not set input format");
            return conv;
        }

    } // namespace detail


    inline void readMoleculesOpenBabel(const std::string &filename,
                              std::function<bool (OpenBabel::OBMol&)> callback)
    {
        OpenBabel::OBMol mol;
        auto conv = detail::makeOBConversion(filename);
        auto notAtEnd = conv.ReadFile(&mol, filename);
        if (!notAtEnd)
            throw std::runtime_error("Could not read molecule");
        while (notAtEnd) {
            if (!callback(mol))
                break;
            notAtEnd = conv.Read(&mol);
        }
    }

    inline OpenBabel::OBMol readSmilesOpenBabel(std::string_view smiles)
    {
        OpenBabel::OBConversion conv;
        conv.SetInFormat("smi");
        OpenBabel::OBMol mol;
        conv.ReadString(&mol, smiles.data());
        return mol;
    }

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
        return get_bonds(mol, atom) | std::views::transform([&mol, atom] (auto bond) { return Kitimar::get_nbr(mol, bond, atom); });
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

    // AromaticLayer

    inline auto is_aromatic_atom(const OpenBabel::OBMol &mol, OpenBabel::OBAtom *atom)
    {
        return atom->IsAromatic();
    }

    inline auto is_aromatic_bond(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        return bond->IsAromatic();
    }



    inline auto is_cyclic_bond(const OpenBabel::OBMol &mol, OpenBabel::OBBond *bond)
    {
        return bond->IsInRing();
    }

} // namespace OpenBabel
