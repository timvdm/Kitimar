#include "Validate.hpp"

void BindingDB_substructure_part_7(OpenBabel::OBMol &mol)
{
    // SMARTS 301 - 350
    //validate_contains<"C[C@]12C[C@H](O)[C@H]3[C@@H](CC[C@H]4C[C@H](O)CC[C@]34C)[C@@H]1CC[C@]2(O)C(=O)CO">(mol); // FIXME: stereo
    validate_contains<"Cc1cc(C=C2C(=O)Nc3ccc(F)cc23)c2ccccc(OCCN3CCCC3C(C)(C)O)c12">(mol);
    validate_contains<"Cc1ccc2[nH]c(=O)c(CC(O)=O)c(-c3ccccc3)c2c1">(mol);
    validate_contains<"Cc1cnc2nc(N)nc(N)c2n1">(mol);
    validate_contains<"ClC1=C(N=C(NC2=CC=CC=C2)S1)C1=NC=CC=C1">(mol);
    validate_contains<"ClC1=C(N=C(NC2=CC=CC=N2)S1)C1=NC=CC=C1">(mol);
    validate_contains<"ClC1=CC=C2NC=C(C2=C1)C1=CCCCC1">(mol);
    validate_contains<"ClC1=CC=C2OC=C(C3=NC4=C(N3)C=CC=C4)C(=O)C2=C1">(mol);
    validate_contains<"ClC1=NNC(=C1)C1=NC2=C(NC1=O)C=CC=C2">(mol);
    validate_contains<"Clc1cccc(Nc2nc(OCC3CCCCC3)c3[nH]cnc3n2)c1">(mol);
    validate_contains<"FC(F)(F)C1=C(Cl)C=CC(NC(=O)NC2=CC=CC=C2)=C1">(mol);
    validate_contains<"FC(F)(F)C1=CC(=CC=C1)N1CCNCC1">(mol);
    validate_contains<"FC(F)(F)C1=CC(NC(=O)NC2=CC=CC=C2)=CC=C1Cl">(mol);
    validate_contains<"FC(F)(F)C1=CN=CC=C1">(mol);
    validate_contains<"FC(F)(F)C1=CN=CC=C1Cl">(mol);
    validate_contains<"FC1=CC=C(CC2=NN=C3N=CC(=NN23)C2=CC=CC=C2)C=C1">(mol);
    validate_contains<"FC1=CC=C(CN2C(NC3CCNCC3)=NC3=C2C=CC=C3)C=C1">(mol);
    validate_contains<"FC1=CC=CC=C1C1=NC2=CC=CC=C2N1">(mol);
    validate_contains<"N(C1=CC=CC=C1)C1=NC=CC(=N1)C1=CN=CS1">(mol);
    validate_contains<"N(C1=CC=CC=C1)C1=NC=NC(=C1)C1=CC=CC=C1">(mol);
    validate_contains<"N(C1=NC=C(O1)C1=CC=CC(=C1)C1=NC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"N(C1=NC=CC=C1)C1=NC(=CC2=NC=CN12)C1=CC=CC=C1">(mol);
    validate_contains<"N1C(=NC(=C1C1=CC=NC=C1)C1=CC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"N1C2=C(C=CC=C2)C2=C1N=C1C=CC=CC1=N2">(mol);
    validate_contains<"N1C2=CC3=C(C=NC=C3)C=C2C2=C1C=CC=C2">(mol);
    validate_contains<"N1C=NC2=C1C1=C(C=CC=C1)N=C2">(mol);
    validate_contains<"N1C=NC2=CC3=NC(=CN=C3C=C12)C1=CC=CC=C1">(mol);
    validate_contains<"N1N=CC2=C(C=CC=C12)C1=NC2=CC=CC=C2N1">(mol);
    validate_contains<"N1N=CC=C1C1=NC2=C(C=CC=C2)N2N=NN=C12">(mol);
    //validate_contains<"N=CC=CC=C1/C(=O)NC(=O)C2=C1C=CC=C2">(mol); // FIXME: stereo
    validate_contains<"NC(=N)C1=CC=C(CNC2=CC=CC=N2)C=C1">(mol);
    validate_contains<"NC(=N)N1CCCC1">(mol);
    validate_contains<"NC(=N)NC1=NC=CC=C1C(N)=O">(mol);
    validate_contains<"NC(=N)NN=C">(mol);
    validate_contains<"NC(=O)C1=CC=CN=C1N">(mol);
    validate_contains<"NC(=O)C1=CC=CN=C1NC1=NC(=CC2=NC=CN12)C1=CC=CC=C1">(mol);
    validate_contains<"NC(=O)C1=CC=CN=C1NC1=NC=CC2=NC=CN12">(mol);
    validate_contains<"NC(=O)C1=CC=CN=C1NC1=NC=CCN1">(mol);
    validate_contains<"NC(=O)C1=CN=C(C=C1)C1=CNC=C1">(mol);
    validate_contains<"NC(=O)C1=CN=C(C=C1)C1=CSC(N)=N1">(mol);
    validate_contains<"NC(=O)C1=CN=C(C=C1)C1=CSC=C1">(mol);
    validate_contains<"NC(=O)C1=CN=C(C=C1)C1=CSC=N1">(mol);
    validate_contains<"NC(=O)C1=NNC2=C1CCC1=C2C=NN1">(mol);
    validate_contains<"NC(=O)NC1=CC=C(OC2=CC(=NC=C2)C(=O)NCCCC(O)=O)C=C1">(mol);
    validate_contains<"NC(=O)NN=C1CCCC2=CC=CN=C12">(mol);
    validate_contains<"NC(=O)NN=C1CCCCC1">(mol);
    //validate_contains<"NC(=S)NN=C1/C(=O)NC2=C1C=CC=C2">(mol); // FIXME: stereo
    validate_contains<"NC(N)N1CCNCC1">(mol);
    validate_contains<"NC1=C(SC(NC2=CC=C(C=C2)N2CCNCC2)=N1)C(=O)C1=C(F)C=CC=C1F">(mol);
    validate_contains<"NC1=C(SC(NC2=CC=CC=C2)=N1)C(=O)C1=CC=CC=C1">(mol);
}
