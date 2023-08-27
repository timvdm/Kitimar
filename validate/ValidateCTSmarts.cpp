#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include "../test/TestData.hpp"
#include <catch2/catch_all.hpp>

using namespace Kitimar;

void Agrafiotis_ABCD_part_1(OpenBabel::OBMol &mol);
void Hicks_and_Jochum_part_1(OpenBabel::OBMol &mol);
void RDKit_smarts_part_1(OpenBabel::OBMol &mol);
void RDKit_smarts_part_2(OpenBabel::OBMol &mol);
void RDKit_smarts_part_3(OpenBabel::OBMol &mol);
void RDKit_smarts_part_4(OpenBabel::OBMol &mol);
void RDKit_smarts_part_5(OpenBabel::OBMol &mol);
void RDKit_smarts_part_6(OpenBabel::OBMol &mol);
void RDKit_smarts_part_7(OpenBabel::OBMol &mol);
void RDKit_smarts_part_8(OpenBabel::OBMol &mol);
void RDKit_smarts_part_9(OpenBabel::OBMol &mol);
void RDMACCS_part_1(OpenBabel::OBMol &mol);
void RDMACCS_part_2(OpenBabel::OBMol &mol);
void RDMACCS_part_3(OpenBabel::OBMol &mol);
void RDMACCS_part_4(OpenBabel::OBMol &mol);
void Rarey_smarts_part_1(OpenBabel::OBMol &mol);
void Rarey_smarts_part_2(OpenBabel::OBMol &mol);
void Rarey_smarts_part_3(OpenBabel::OBMol &mol);
void Rarey_smarts_part_4(OpenBabel::OBMol &mol);
void Rarey_smarts_part_5(OpenBabel::OBMol &mol);
void Rarey_smarts_part_6(OpenBabel::OBMol &mol);
void Rarey_smarts_part_7(OpenBabel::OBMol &mol);
void Rarey_smarts_part_8(OpenBabel::OBMol &mol);
void Rarey_smarts_part_9(OpenBabel::OBMol &mol);
void Rarey_smarts_part_10(OpenBabel::OBMol &mol);
void Rarey_smarts_part_11(OpenBabel::OBMol &mol);
void Rarey_smarts_part_12(OpenBabel::OBMol &mol);
void Rarey_smarts_part_13(OpenBabel::OBMol &mol);
void Rarey_smarts_part_14(OpenBabel::OBMol &mol);
void Rarey_smarts_part_15(OpenBabel::OBMol &mol);
void Rarey_smarts_part_16(OpenBabel::OBMol &mol);
void Rarey_smarts_part_17(OpenBabel::OBMol &mol);
void Rarey_smarts_part_18(OpenBabel::OBMol &mol);
void Rarey_smarts_part_19(OpenBabel::OBMol &mol);
void Rarey_smarts_part_20(OpenBabel::OBMol &mol);
void Rarey_smarts_part_21(OpenBabel::OBMol &mol);
void Rarey_smarts_part_22(OpenBabel::OBMol &mol);
void Rarey_smarts_part_23(OpenBabel::OBMol &mol);
void Rarey_smarts_part_24(OpenBabel::OBMol &mol);
void Rarey_smarts_part_25(OpenBabel::OBMol &mol);

TEST_CASE("Substructure Quaery Collection")
{
    OpenBabelSmilesMolSource source{chembl_smi_filename("100K")};
    auto i = 0;
    for (auto mol : source.molecules()) {
        std::cout << "Molecule: " << ++i << " -- " << writeSmiles(mol) << '\n';
        Agrafiotis_ABCD_part_1(mol);
        Hicks_and_Jochum_part_1(mol);
        RDKit_smarts_part_1(mol);
        RDKit_smarts_part_2(mol);
        RDKit_smarts_part_3(mol);
        RDKit_smarts_part_4(mol);
        RDKit_smarts_part_5(mol);
        RDKit_smarts_part_6(mol);
        RDKit_smarts_part_7(mol);
        RDKit_smarts_part_8(mol);
        RDKit_smarts_part_9(mol);
        RDMACCS_part_1(mol);
        RDMACCS_part_2(mol);
        RDMACCS_part_3(mol);
        RDMACCS_part_4(mol);
        Rarey_smarts_part_1(mol);
        Rarey_smarts_part_2(mol);
        Rarey_smarts_part_3(mol);
        Rarey_smarts_part_4(mol);
        Rarey_smarts_part_5(mol);
        Rarey_smarts_part_6(mol);
        Rarey_smarts_part_7(mol);
        Rarey_smarts_part_8(mol);
        Rarey_smarts_part_9(mol);
        Rarey_smarts_part_10(mol);
        Rarey_smarts_part_11(mol);
        Rarey_smarts_part_12(mol);
        Rarey_smarts_part_13(mol);
        Rarey_smarts_part_14(mol);
        Rarey_smarts_part_15(mol);
        Rarey_smarts_part_16(mol);
        Rarey_smarts_part_17(mol);
        Rarey_smarts_part_18(mol);
        Rarey_smarts_part_19(mol);
        Rarey_smarts_part_20(mol);
        Rarey_smarts_part_21(mol);
        Rarey_smarts_part_22(mol);
        Rarey_smarts_part_23(mol);
        Rarey_smarts_part_24(mol);
        Rarey_smarts_part_25(mol);
    }
}
