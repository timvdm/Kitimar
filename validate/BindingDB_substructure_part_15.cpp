#include "Validate.hpp"

void BindingDB_substructure_part_15(OpenBabel::OBMol &mol)
{
    // SMARTS 701 - 750
    validate_contains<"C1CCC(CC1)C2CCCC2C3=CC=CC=C3">(mol);
    validate_contains<"C1CCC(CC1)C2CCNCC2">(mol);
    validate_contains<"C1CCC(CC1)NC2=C3N=CNC3=NC=N2">(mol);
    validate_contains<"C1CCC2(C1)CCCC3=C2C4=C(C=CC=C4)C=C3">(mol);
    validate_contains<"C1CCC2(CC1)CCC3=CC4=C(C=CC=C4)C=C3C2">(mol);
    validate_contains<"C1CCC2(CC1)CCCC2">(mol);
    validate_contains<"C1CCC2=C(C1)C(C3=CC=CC=C3)C4=C(CCCC4)N2C5=CC=CC=C5">(mol);
    validate_contains<"C1CCC2=C(C1)NC3=C2C=CC=C3">(mol);
    validate_contains<"C1CCC2=CC3=C(C=CC=C3)C=C2C1">(mol);
    validate_contains<"C1CCC2=CC3=C(CC4=CC5=C(C=CC=C5)C=C34)C=C2C1">(mol);
    validate_contains<"C1CCC2=CC3=C(OCCO3)C=C2C1">(mol);
    validate_contains<"C1CCC2=NC3=CC=CC=C3N2C1">(mol);
    validate_contains<"C1CCC2C(C1)C(C3CCCCC3N2C4=CC=CC=C4)C5=CC=CC=C5">(mol);
    validate_contains<"C1CCC2C(C1)CCC3=CC=CC=C23">(mol);
    validate_contains<"C1CCC2C(C1)CCC3C2CCC4C5CCCCC5CCC34">(mol);
    validate_contains<"C1CCC2CCCCCC2C1">(mol);
    validate_contains<"C1CCC2OCCC2C1">(mol);
    validate_contains<"C1CCC2OCCCC2C1">(mol);
    validate_contains<"C1CCC2SCCCC2C1">(mol);
    validate_contains<"C1CCC=CC1">(mol);
    //validate_contains<"C1CCCC1.N2C=CC=C2.C3=CC4=C(C=C3)C=CC=C4">(mol); // FIXME: components
    //validate_contains<"C1CCCCC1.C2=CC3=C(C=C2)C=CC=C3">(mol); // FIXME: components
    validate_contains<"C1CCN(CC1)C2(CCCCC2)C3=CC=CC=C3">(mol);
    validate_contains<"C1CCN(CC1)C2=CC=CC=C2">(mol);
    validate_contains<"C1CCN(CC1)C2=CC=NC=C2">(mol);
    validate_contains<"C1CCN(CC1)C2=NC=NC3=C2C=CN3">(mol);
    validate_contains<"C1CCN2CCC3=CC=CC=C3C2C1">(mol);
    validate_contains<"C1CN(CCN1)C2=CC=CC3=C2C=CN=C3">(mol);
    validate_contains<"C1CN(CCN1C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"C1CN2C(CC1C3=CC=CC=C3)=NC4=C2C=CC=C4">(mol);
    validate_contains<"C1CN2CC3=C(CNCC3)N=C2C1">(mol);
    validate_contains<"C1CN2CC=CN=C2C1">(mol);
    validate_contains<"C1CN2CCC1CC2">(mol);
    //validate_contains<"C1CN2CC[C@H]1CC2">(mol); // FIXME: stereo
    validate_contains<"C1CN2N=NN=C2C3=CC=CC=C3N1">(mol);
    validate_contains<"C1CN=C(N1)C2=CC3=C(O2)C=CC=C3">(mol);
    validate_contains<"C1CN=C(N1)C2=CC=C(C=C2)C3=CC=C(O3)C4=CC=CC=C4">(mol);
    validate_contains<"C1CN=C(N1)C2=CC=C(C=C2)C3=CC=CO3">(mol);
    validate_contains<"C1CN=C(NN=CC2=CC=CC=C2)N1">(mol);
    validate_contains<"C1CN=CN1">(mol);
    validate_contains<"C1CNC(C1)C2CNCC2C3CCCN3">(mol);
    validate_contains<"C1CNC2=CC=CC=C2SC1">(mol);
    validate_contains<"C1CNC2=CC=CN=C2C1">(mol);
    validate_contains<"C1CNC2CCCNC2C1">(mol);
    validate_contains<"C1CNc2ccccc2CN1">(mol);
    validate_contains<"C1COC2=CC(=CC=C2O1)C3=C(NC=N3)C4=NC=CC=C4">(mol);
    validate_contains<"C1COCO1">(mol);
    validate_contains<"C1N=CC=[S]1">(mol);
    validate_contains<"C1NC(=C2CNc3ccccc23)c4ccccc14">(mol);
    validate_contains<"C1NC2=C(C=CC=C2)C1C3=CC=CC=C3">(mol);
}
