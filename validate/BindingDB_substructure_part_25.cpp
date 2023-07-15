#include "Validate.hpp"

void BindingDB_substructure_part_25(OpenBabel::OBMol &mol)
{
    // SMARTS 1201 - 1250
    validate_contains<"CC1CCC2=CC=CC=C2C1">(mol);
    validate_contains<"CC1CCC2C1C(O)CC3C2CCC4=CC(O)=CC=C34">(mol);
    validate_contains<"CC1CCC2C3CCC4=C(C=CC(O)=C4)C3CCC12O">(mol);
    validate_contains<"CC1CCC2CCCCC2C1">(mol);
    validate_contains<"CC1CCC=CC(C)=CC=CC=C1">(mol);
    validate_contains<"CC1CCCC(=O)N1">(mol);
    validate_contains<"CC1CCCC(C)N1C">(mol);
    //validate_contains<"CC1CCCC1.CC2CCCC2">(mol); // FIXME: components
    validate_contains<"CC1CCCC1CCC2CCCC3CCCCC23">(mol);
    validate_contains<"CC1CCCC2CCCCC12">(mol);
    //validate_contains<"CC1CCCCC1.C2CCC3(CC2)CCC4(CCCCC4)CC3">(mol); // FIXME: components
    validate_contains<"CC1CCCCC1OC2=C(O)C(=CC=C2)C3=CC4=CC(C(N)=N)=C(F)C=C4N3">(mol);
    validate_contains<"CC1CCCCCC(O)C(=O)C(O)CC(=O)O1">(mol);
    validate_contains<"CC1CCCCN1">(mol);
    validate_contains<"CC1CCN=C1">(mol);
    //validate_contains<"CC1CCNC1C=C2/C(=O)CC3=CC=CC=C23">(mol); // FIXME: stereo
    validate_contains<"CC1CN(C)CC1C">(mol);
    validate_contains<"CC1CN(C)CC1C2=CN=CC=N2">(mol);
    validate_contains<"CC1CN(CC(C)N1)C2=CC3=C(C=C2F)C(=O)C(=CN3C(C)(C)C)C(O)=O">(mol);
    validate_contains<"CC1CNC(CN1)C=O">(mol);
    validate_contains<"CC1CNCC(N1)C=O">(mol);
    validate_contains<"CC1CNCC1CN">(mol);
    validate_contains<"CC1NC(=O)NC1=O">(mol);
    validate_contains<"CC1NC(=O)NC=C1">(mol);
    validate_contains<"CC1NC2=CC=CC=C2OC1O">(mol);
    validate_contains<"CC1NCC(N)C1O">(mol);
    validate_contains<"CC1NNNN1">(mol);
    validate_contains<"CC1OC(C(C)C1C)N2C=NC3=C2N=CNCC3O">(mol);
    validate_contains<"CC1OC2=CC=CC=C2N(C)C1=O">(mol);
    validate_contains<"CC1OC2C(C2C1(N)C(O)=O)C(O)=O">(mol);
    //validate_contains<"CC=C/C=CC.C=CC=C/C=C.C1=CC=CC=C1">(mol); // FIXME: components
    //validate_contains<"CC=C1/C(=O)NC(=O)C2=C1C=CC=C2">(mol); // FIXME: stereo
    //validate_contains<"CC=C1/C2CC3=C(C=CC(=O)N3)C1(N)CC(C)=C2">(mol); // FIXME: stereo
    validate_contains<"CC=C1C(=O)NC2=C1C=CC=C2">(mol);
    validate_contains<"CC=CCC=CCC1CC2=CC=CC=C2C1C">(mol);
    validate_contains<"CCC(=N)NC(=S)NC">(mol);
    validate_contains<"CCC(=O)N1CCCC1C(=O)NC">(mol);
    validate_contains<"CCC(=O)NC(=S)NC">(mol);
    //validate_contains<"CCC(=O)N[C@@H](CC(O)=O)C(=O)NC">(mol); // FIXME: stereo
    validate_contains<"CCC(=O)OCCC(C)C">(mol);
    validate_contains<"CCC(C)(C)C(=O)C(=O)N1CCCC1C(=O)OCCCC2=CC=CC=C2">(mol);
    validate_contains<"CCC(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12">(mol);
    validate_contains<"CCC(C)C(C(=O)N1CCC(F)C1)N">(mol);
    validate_contains<"CCC(C)C1=CC=CC(C)=C1">(mol);
    validate_contains<"CCC(C1=CC(=O)NC(SC(C)C)=N1)C2=CC=CC3=C2C=CC=C3">(mol);
    validate_contains<"CCC(CC)C1C=C(C)CC(NC(=O)CCOCCOCCOCCOCCN)C1NC(C)=O">(mol);
    validate_contains<"CCC(CC)NC1=NC(C)=NC2=C(C(C)=NN12)C3=C(C)N=C(C=C3C)N(C)C">(mol);
    validate_contains<"CCC(CC)OC1=C(C)C(CC2=C(C)C=C(C)C=C2C)=NC(C)=C1">(mol);
    validate_contains<"CCC(Cc1ccc(OC)c(CNC(=O)c2ccc(cc2F)C(F)(F)F)c1)C(O)=O">(mol);
    validate_contains<"CCC(N)=N">(mol);
}
