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

    inline auto chembl_serialized_filename(CTLayout::StructMolecule)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMolecule";
    }

    inline auto chembl_serialized_filename(CTLayout::StructMoleculeIncident)
    {
        return std::string(KITIMAR_DATA_DIR) + "/chembl_32.StructMoleculeIncident";
    }



    /*
    template<typename FromLayout, typename ToLayout>
    void serialize_chembl()
    {
        std::vector<std::byte> idata, odata;
        std::ifstream ifs(chembl_serialized_filename(FromLayout{}), std::ios_base::binary | std::ios_base::in);
        std::ofstream ofs(chembl_serialized_filename(ToLayout{}), std::ios_base::binary | std::ios_base::out);
        assert(ofs);
        while (ifs) {
            auto [ok, mol] = deserialize<FromLayout>(ifs, idata);
            if (!ok)
                break;
            serialize<ToLayout>(mol, ofs, odata);
        }

        ofs.write(reinterpret_cast<const char*>(&LayoutEnd), LayoutSize::size());
    }
    */

} // namespace Kitimar
