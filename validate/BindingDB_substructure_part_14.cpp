#include "Validate.hpp"

void BindingDB_substructure_part_14(OpenBabel::OBMol &mol)
{
    // SMARTS 651 - 700
    validate_contains<"C1C=C1">(mol);
    validate_contains<"C1C=CC2=C(C=CC=C2)C1NC3=CC=CC=C3">(mol);
    validate_contains<"C1C=CC2=C1C=CC=C2">(mol);
    validate_contains<"C1C=CC2=C3C=CC=CC3=CC=C12">(mol);
    validate_contains<"C1C=CC2=CC3=CC4=C5C(C=CC4)=CC6=CC7=CC=CC8=C7C9=CC(C=C1C2C=C3C5=C69)=C8">(mol);
    validate_contains<"C1C=NC(=C1)C2=C(CNC2)C3=CC=CN3">(mol);
    validate_contains<"C1C=NC(N=NC2=CC=CC=C2)=[S]1">(mol);
    validate_contains<"C1C=NOC12CCNCC2">(mol);
    validate_contains<"C1CC(C2=CC=CC=C2)C3=CC=CC=C3C1">(mol);
    validate_contains<"C1CC(C=N1)C2=NC=NC=C2">(mol);
    validate_contains<"C1CC(CCN1)(c2ccccc2)c3ccc(cc3)-c4cn[nH]c4">(mol);
    //validate_contains<"C1CC(CCN1)=C(/C2=CC=CC=C2)C3=CC=CC=C3">(mol); // FIXME: stereo
    validate_contains<"C1CC(CCN1)C2=CC=CC=C2">(mol);
    validate_contains<"C1CC(CCN1)C2=CC=CN=N2">(mol);
    validate_contains<"C1CC(CCN1)C2=COC=N2">(mol);
    validate_contains<"C1CC(CCN1)C2=NNC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"C1CC(CCN1)CC2=CC(OCC3=CC=CC=C3)=CC(=C2)C4=CC=CC=C4">(mol);
    validate_contains<"C1CC(CCS1)C2=CC=CC=C2">(mol);
    validate_contains<"C1CC(CN1)C2=CC=CC=C2">(mol);
    validate_contains<"C1CC2(C(N1)NC3=C2C=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"C1CC2=C(C1)C3=C(C=CC=C3)C=C2">(mol);
    validate_contains<"C1CC2=C(C1)C3=C(N2)NC4=C3CCCC4">(mol);
    validate_contains<"C1CC2=C(C1)C=C3C=CC=CC3=C2">(mol);
    validate_contains<"C1CC2=C(CN1)C=CS2">(mol);
    validate_contains<"C1CC2=C(N1)C3=C(C4=C2CNC4)C5=CC=CC=C5N3">(mol);
    validate_contains<"C1CC2=C(N1)C3=C(CCN3)C4=C2CNC4">(mol);
    validate_contains<"C1CC2=C(NC3=CC=CC=C23)C(N1)C4=CC=CC=C4">(mol);
    validate_contains<"C1CC2=CC3=C(C=CN3)C=C2C1">(mol);
    validate_contains<"C1CC2=CC3=C(CCO3)C=C2O1">(mol);
    validate_contains<"C1CC2=CC3=C(OCC3)C=C2O1">(mol);
    validate_contains<"C1CC2=CC3=CC=C4C5=CC=CC=C5C=CC4=C3C=C2C1">(mol);
    validate_contains<"C1CC2=CC=CC3=C2C(C1)=CC=C3">(mol);
    validate_contains<"C1CC2=CC=CC3=C2C1=CC=C3">(mol);
    validate_contains<"C1CC2=CC=CC=C2C(O1)C3=C4C=CC=CC4=CN3">(mol);
    validate_contains<"C1CC2=CC=CC=C2C3=CC=CN13">(mol);
    validate_contains<"C1CC2=CC=CC=C2CCN1">(mol);
    validate_contains<"C1CC2=NN=CN2C3=C(C1)C=CC=C3">(mol);
    validate_contains<"C1CC2C(C1)C3=CC4=C(C=CC=C4)C=C23">(mol);
    //validate_contains<"C1CC2CC3=C(CC2C1)C=C4C=CC=CC4=C3.C(CCC5CCC6C7CCCC7CCC6C5)CC8CCCC8">(mol); // FIXME: components
    validate_contains<"C1CC2CC3=C(CC2C1)C=CN3">(mol);
    validate_contains<"C1CC2CCC1C2">(mol);
    validate_contains<"C1CC2CCC1N2">(mol);
    validate_contains<"C1CC2CCC3=C(C=CC4=CC=CC=C34)C2C1">(mol);
    validate_contains<"C1CC2CCC3C(CCC4=C3C=CC=C4)C2C1">(mol);
    validate_contains<"C1CC2CCC3CCCC4CCC(C1)C2C34">(mol);
    validate_contains<"C1CC2CCCC2C1">(mol);
    validate_contains<"C1CC2CCCC3CCCC(C1)C23">(mol);
    validate_contains<"C1CC2CCNC2N1">(mol);
    validate_contains<"C1CC2CCOC2N1">(mol);
    validate_contains<"C1CCC(C1)CCC2=CC=C3C=CC=CC3=C2">(mol);
}
