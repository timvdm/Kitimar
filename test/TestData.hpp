#pragma once

#include <Kitimar/Util/Util.hpp>
#include <Kitimar/CTLayout/Molecule.hpp>

namespace Kitimar {

    namespace detail {

        inline auto addDot(const std::string &s)
        {
            if (s.empty() || s.starts_with('.'))
                return s;
            return "." + s;
        }

    }

    inline auto chembl_smi_filename(const std::string &suffix = {})
    {
        return fmt::format("{}/chembl_32{}.smi", KITIMAR_DATA_DIR, detail::addDot(suffix));
    }

    inline auto chembl_rdkit_filename(const std::string &suffix = {})
    {
        return fmt::format("{}/chembl_32{}.pkl", KITIMAR_DATA_DIR, detail::addDot(suffix));
    }



    inline auto chembl_serialized_filename(CTLayout::Vector<CTLayout::StructMolecule>, const std::string &suffix = {})
    {
        return fmt::format("{}/chembl_32{}.StructMolecule", KITIMAR_DATA_DIR, detail::addDot(suffix));
    }

    inline auto chembl_serialized_filename(CTLayout::Vector<CTLayout::StructMoleculeIncident>, const std::string &suffix = {})
    {
        return fmt::format("{}/chembl_32{}.StructMoleculeIncident", KITIMAR_DATA_DIR, detail::addDot(suffix));
    }

    inline auto chembl_serialized_filename(CTLayout::Vector<CTLayout::ListMoleculeIncident>, const std::string &suffix = {})
    {
        return fmt::format("{}/chembl_32{}.ListMoleculeIncident", KITIMAR_DATA_DIR, detail::addDot(suffix));
    }

    inline auto chembl_serialized_filename(CTLayout::TypeMolecules, const std::string &suffix = {})
    {
        return fmt::format("{}/chembl_32{}.TypeMolecule", KITIMAR_DATA_DIR, detail::addDot(suffix));
    }


} // namespace Kitimar
