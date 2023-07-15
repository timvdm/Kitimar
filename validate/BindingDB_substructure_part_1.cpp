#include "Validate.hpp"

void BindingDB_substructure_part_1(OpenBabel::OBMol &mol)
{
    // SMARTS 1 - 50
    validate_contains<"BrC1=C2NCCCCNC3=CC=CC(NC(N=C1)=N2)=C3">(mol);
    validate_contains<"C(CC1=CC=CC=C1)CC1=C(C=CC=C1)C1CCCCC1">(mol);
    validate_contains<"C(N1CCC(CC1)NC1=CC=NC2=C1C=CC=C2)C1=CC=CC=C1">(mol);
    validate_contains<"C(NC1=CC=CC2=C1CCC2)C1=CC=CC=C1">(mol);
    validate_contains<"C(OC1=CC2=C(NCC2)C=C1)C1=NC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"C(c1ccncc1)c1cc(-c2ccccc2)c2ncccc2c1">(mol);
    //validate_contains<"C.CN1CC[C@]23[C@H]4OC5=C(O)C=CC(C[C@@H]1[C@]2(O)CCC4=O)=C35">(mol); // FIXME: components
    validate_contains<"C1=CC2=C(C=C1)C1=C(C=C2)N=CC(=C1)C1=CC=CN=C1">(mol);
    validate_contains<"C1=CC2=C(C=C1)N=CC=N2">(mol);
    validate_contains<"C1=CC=C(C=C1)C1=CC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"C1=CC=C(C=C1)C1=CC=CC2=C1C=CC=C2">(mol);
    //validate_contains<"C1=CC=CC=C1.C1=CC2=C(C=C1)C=CC=C2">(mol); // FIXME: components
    validate_contains<"C1=CN2C=CC=NC2=C1">(mol);
    validate_contains<"C1=NC2C=CC=NC2=C1">(mol);
    validate_contains<"C1=NC2C=NC=CC2=C1">(mol);
    validate_contains<"C1=NC2C=NC=NC2=C1">(mol);
    validate_contains<"C1C2=C(OC3=C1C1=C(C=CC=C1)C=C3)N=CN=C2">(mol);
    validate_contains<"C1C2=CC=C(CC3=CC=C(CC4=CC=C(CC5=CC=C1N5)N4)N3)N2">(mol);
    validate_contains<"C1C=NC2C=NC=CN12">(mol);
    validate_contains<"C1CC2**C3=C(C=CC=C3)C2CN1">(mol);
    validate_contains<"C1CCC(=CC1)N1CCOCC1">(mol);
    validate_contains<"C1CCC(C1)C1=CC=CC2=C1N=C(NC1=CC=C(C=C1)N1CCNCC1)N=C2">(mol);
    validate_contains<"C1CCC(CC1)C1CCCCC1">(mol);
    validate_contains<"C1CCC(NC1)C1=CC=CC=C1">(mol);
    validate_contains<"C1CCC2CC3=C(CC2C1)NC=C3">(mol);
    //validate_contains<"C1CCCC1.C1CC2CCCC2C1">(mol); // FIXME: components
    validate_contains<"C1CCCCCCCC1">(mol);
    validate_contains<"C1CCN(C1)C1=CSC=N1">(mol);
    validate_contains<"C1CCNC2=CC=NC(NC3=CC(NC1)=CC=C3)=N2">(mol);
    validate_contains<"C1CN2C(CN1)**C1=C2C=CC=C1">(mol);
    validate_contains<"C1CN2CCC3=C(C=CC=C3)C2CN1">(mol);
    validate_contains<"C1CN2N=CNC2S1">(mol);
    validate_contains<"C1CNC(C1)C1=CC=CC=C1">(mol);
    validate_contains<"C1CNC(C1)C1=CC=CC=N1">(mol);
    validate_contains<"C1CNC(CN1)C1=CC=CC=C1">(mol);
    validate_contains<"C1CNC(CN1)C1=CC=NC=C1">(mol);
    validate_contains<"C1CNC2C=CSC2C1">(mol);
    validate_contains<"C1COC2=NC=CN2C1">(mol);
    validate_contains<"C1COCCO1">(mol);
    validate_contains<"C1C[N+]2=C(NC=N2)S1">(mol);
    validate_contains<"C1NC2=C(NC3=C1C=CN=N3)C=CC=C2">(mol);
    validate_contains<"C1OC2=C(C=CC=C2)C2=C1NC=N2">(mol);
    validate_contains<"C1OC2=C(O1)C(=CC=C2)C1=CC(NC2=CC=CC=C2)=NC=N1">(mol);
    validate_contains<"C1OCC2=COCC12">(mol);
    validate_contains<"C1OCC2COCC12">(mol);
    validate_contains<"C1OCN2CC=CC3C=CCN1C23">(mol);
    validate_contains<"C=C1NC(=O)SC1CC1=CC=C(OCOC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"C=CC1=NC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"CC(=C)C=C">(mol);
    //validate_contains<"CC(=C/C1=CC=C(C=C1)C(O)=O)C1=CC2=C(C=C1)C(C)(C)CCC2(C)C">(mol); // FIXME: stereo
}
