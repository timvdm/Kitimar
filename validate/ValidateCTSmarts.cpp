#include "Validate.hpp"
#include "../test/TestData.hpp"
#include <catch2/catch_all.hpp>

using namespace Kitimar;

//
// Basic
//

#ifdef KITIMAR_WITH_VALIDATION_BASIC

template<Molecule::Molecule Mol>
void validate_basic(Mol &mol);

template<Toolkit::Id ToolkitId>
void toolkit_validate_basic()
{
    auto source = Toolkit::smilesMolSource<ToolkitId>(chembl_smi_filename("100K"));
    auto i = 0;
    for (auto mol : source.molecules()) {
        std::cout << "Molecule: " << ++i << " -- " << Toolkit::writeSmiles<ToolkitId>(*mol) << '\n';
        validate_basic(*mol);
    }
}

TEST_CASE("rdkit_basic")
{
    toolkit_validate_basic<Toolkit::rdkit>();
}


TEST_CASE("openbabel_basic")
{
    toolkit_validate_basic<Toolkit::openbabel>();
}

#endif // KITIMAR_WITH_VALIDATION_BASIC

//
// Substructure Query Collection
//

#ifdef KITIMAR_WITH_VALIDATION_SQC

template<Molecule::Molecule Mol>
void validate_parts(Mol &mol);

template<Toolkit::Id ToolkitId>
void toolkit_validate_parts()
{
    auto source = Toolkit::smilesMolSource<ToolkitId>(chembl_smi_filename());
    auto i = 0;
    for (auto mol : source.molecules()) {
        std::cout << "Molecule: " << ++i << " -- " << Toolkit::writeSmiles<ToolkitId>(*mol) << '\n';
        validate_parts(*mol);
    }
}

TEST_CASE("rdkit_sqc")
{
    toolkit_validate_parts<Toolkit::rdkit>();
}

TEST_CASE("openbabel_sqc")
{
    toolkit_validate_parts<Toolkit::openbabel>();
}

#endif // KITIMAR_WITH_VALIDATION_SQC
