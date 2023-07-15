#include "Validate.hpp"

void BindingDB_substructure_part_46(OpenBabel::OBMol &mol)
{
    // SMARTS 2251 - 2300
    validate_contains<"O=C1CC2=CC=CC=C2C(=O)N1">(mol);
    validate_contains<"O=C1CC2CCC1C(N2C3=CC=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"O=C1CC2SCC=CC12">(mol);
    validate_contains<"O=C1CC=CC(=O)N1">(mol);
    validate_contains<"O=C1CCC(=O)C2=C1NC=C2">(mol);
    validate_contains<"O=C1CCC(=O)N1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CCC2=CC=CC=C12">(mol);
    validate_contains<"O=C1CCC2=CC=CC=C2C1">(mol);
    validate_contains<"O=C1CCC2CC3=CC4=C(OCO4)C=C3CC12">(mol);
    validate_contains<"O=C1CCC2CCCCC2C1">(mol);
    validate_contains<"O=C1CCCC(=O)C1">(mol);
    validate_contains<"O=C1CCCC(=O)C1C(=O)C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CCCC1CC2CCCCC2">(mol);
    validate_contains<"O=C1CCCC2=C1C(C3=CC=CC=C3)C4=C(C2)C=CC=C4">(mol);
    validate_contains<"O=C1CCCC2=C1C=CN2">(mol);
    validate_contains<"O=C1CCCC2CCCN12">(mol);
    validate_contains<"O=C1CCCC2OCCN12">(mol);
    validate_contains<"O=C1CCCCC1">(mol);
    validate_contains<"O=C1CCCCCN1">(mol);
    validate_contains<"O=C1CCCCN1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CCCCO1">(mol);
    validate_contains<"O=C1CCCN1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CCN1">(mol);
    validate_contains<"O=C1CCO1">(mol);
    validate_contains<"O=C1CN(CC2N1CCC3=C2C=CC=C3)C(=O)C4CCCCC4">(mol);
    validate_contains<"O=C1CN=C(C2=CC=CC=C2)C3=CC=CC=C3N1">(mol);
    validate_contains<"O=C1CNC(=O)CN1">(mol);
    validate_contains<"O=C1COCC(=O)N1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1COCCN1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CSC(=S)N1">(mol);
    validate_contains<"O=C1N(C2=CC=C(C=C2C1=O)S(=O)(=O)N3CCCC3COC4=CC=CC=C4)C5=CC=CC=C5">(mol);
    validate_contains<"O=C1N(C2CCC(=O)NC2=O)C(=O)C3=CC=CC=C13">(mol);
    validate_contains<"O=C1N(CC2=CC=CC=C2)C3=CC=C(C=C3C1=O)S(=O)(=O)N4CCCC4">(mol);
    validate_contains<"O=C1N(CC2=CC=CC=C2)C3=CC=C(C=C3C1=O)S(=O)(=O)N4CCCC4COC5=CC=CC=C5">(mol);
    validate_contains<"O=C1N([NH+]=CC2=CC=CC=C2)C(=NC3=CC=CC=C13)C4=CC=CC=C4">(mol);
    //validate_contains<"O=C1N=C2C=CC=CC2=CC1=C3/NC4=CC=CC=C4N3">(mol); // FIXME: stereo
    validate_contains<"O=C1N=CN=C2NNC=C12">(mol);
    validate_contains<"O=C1NC(=NC2=C1C=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"O=C1NC(=O)C(C1NC2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"O=C1NC(=O)C(C=C1)C2CCCO2">(mol);
    validate_contains<"O=C1NC(=O)C2=C3C=CC=CC3=C4C(=O)NC(=O)C5=C6C=CC=CC6=C1C2=C45">(mol);
    validate_contains<"O=C1NC(=O)C2=CC=CC=C12">(mol);
    validate_contains<"O=C1NC(=O)C2=CC=CC=C2C1=O">(mol);
    validate_contains<"O=C1NC(=O)C2=CC=CC=C2N1">(mol);
    validate_contains<"O=C1NC(=O)N(N=C1)C2=CC=CC=C2">(mol);
    validate_contains<"O=C1NC(=O)N2CCC3=C(NC4=C3C=CC=C4)C12">(mol);
    validate_contains<"O=C1NC2=C(C(=O)C3=CC=CC=C23)C4=CC=CC=C14">(mol);
    validate_contains<"O=C1NC2=C(C3=C(CCCC3)S2)C(=O)N1">(mol);
    validate_contains<"O=C1NC2=C(C3=C(CCCC3)S2)C(=O)N1C4=CC=CC=C4">(mol);
    validate_contains<"O=C1NC2=C(C=C1)C=NC(NC3CCNCC3)=N2">(mol);
}
