#pragma once

#include <Kitimar/CTLayout/Molecule.hpp>

namespace Kitimar {


    inline auto chembl_smi_filename()
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.smi";
    }

    inline auto chembl_smi_filename(const std::string &infix)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32." + infix + ".smi";
    }

    inline auto chembl_rdkit_filename()
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.pkl";
    }

    inline auto chembl_serialized_filename(CTLayout::Vector<CTLayout::StructMolecule>)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMolecule";
    }

    inline auto chembl_serialized_filename(CTLayout::Vector<CTLayout::StructMoleculeIncident>)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMoleculeIncident";
    }

    inline auto chembl_serialized_filename(CTLayout::Vector<CTLayout::AtomTypeMolecule>)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.AtomTypeMolecule";

    }


} // namespace Kitimar
