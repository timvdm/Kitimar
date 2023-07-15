#include "Validate.hpp"

void BindingDB_substructure_part_4(OpenBabel::OBMol &mol)
{
    // SMARTS 151 - 200
    validate_contains<"CC1N2N=CN(C)C2=NC2=C1C(C)C1=C(O2)C=CC2=C1C=CC=C2">(mol);
    validate_contains<"CC1NNC(C)=C1">(mol);
    validate_contains<"CC=CC1=NC2=CC=CC=C2C(=O)N1C">(mol);
    //validate_contains<"CC=CNC=NCNC[C@@H](C)C=C.CN1CCN(CC1)C1=CC=C(NC2=NC(NC[C@H]3C[C@H](CC(N)=O)C=C3)=C(F)C=N2)C=C1C.CN1CCN(CC1)C1=CC=C(NC2=NC(N[C@@H]3[C@H]4C[C@@H](C=C4)C3C(N)=O)=C(F)C=N2)C=C1C.CN1CCN(CC1)C1=CC=C(NC2=NC(N[C@@H]3[C@H]4C[C@@H](C=C4)C3C(N)=O)=C(F)C=N2)C=C1C.CN1CCN(CC1)C1=CC=C(NC2=NC(N[C@@H]3[C@H]4C[C@@H](C=C4)C3C(N)=O)=C(F)C=N2)C=C1C.CN1CCN(CC1)C1=CC=C(NC2=NC(N[C@@H]3[C@H]4C[C@@H](C=C4)C3C(N)=O)=C(F)C=N2)C=C1C.CN1CCN(CC1)C1=CC=C(NC2=NC(N[C@@H]3[C@H]4C[C@@H](C=C4)C3C(N)=O)=C(F)C=N2)C=C1C.CN1CCN(CC1)C1=CC=C(NC2=NC(N[C@@H]3[C@H]4C[C@@H](C=C4)C3C(N)=O)=C(F)C=N2)C=C1C">(mol); // FIXME: components
    validate_contains<"CCC(=O)NC1=NC=CS1">(mol);
    validate_contains<"CCC(C)C(C(CC(=O)N1CCCC1C(OC)C(C)C(=O)NC(CC1=CC=CC=C1)C1=NC=CS1)OC)N(C)C(=O)C(NC(=O)C(C(C)C)N(C)C)C(C)C">(mol);
    validate_contains<"CCC(C)C1=CC=C(C=C1)C(C)C(C)C">(mol);
    //validate_contains<"CCC(C)CC(C)CCCCCCCCC(=O)N[C@@H]1C[C@@H](O)[C@@H](O)NC(=O)[C@@H]2[C@@H](O)CCN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC1=O)[C@@H](C)O)[C@@H](O)[C@H](O)c1ccc(O)c(N)c1)[C@@H](O)CC(N)=O">(mol); // FIXME: stereo
    validate_contains<"CCC(C1=CC(=O)NC(SC(C)C)=N1)C1=CC=CC2=C1C=CC=C2">(mol);
    validate_contains<"CCC(CO)NC1=NC(NCC2=CC=CC=C2)=C2N=CN(C(C)C)C2N1">(mol);
    validate_contains<"CCC(CO)NC1=NC(NCC2=CC=CC=C2)=C2N=CNC2N1">(mol);
    validate_contains<"CCC(N)C1CCCC(C)C1N">(mol);
    validate_contains<"CCC(N)C1CCCCC1N">(mol);
    validate_contains<"CCC(N)CC(C)N">(mol);
    //validate_contains<"CCC.CCC(C)C(C)C">(mol); // FIXME: components
    //validate_contains<"CCC.CCCC.CCCC(C)C">(mol); // FIXME: components
    validate_contains<"CCC1=C(C)N(C)C=C1C1=CC=CC=C1">(mol);
    validate_contains<"CCC1=C(OC2=CC(=C(C=C2)C#N)C(F)(F)F)C=CC=C1">(mol);
    validate_contains<"CCC1=CC=CC(CC)=C1NC(=O)N1CC2=C(C1)C(NC(=O)C1=CC=C(C=C1)N1CCN(C)CC1)=NN2">(mol);
    validate_contains<"CCC1=CN=C(CSC2=CN=C(NC(C)=O)S2)O1">(mol);
    validate_contains<"CCC1=NC2=C(C=CC=C2)C(OC2CCCCC2)=N1">(mol);
    validate_contains<"CCC1=NN(C(NC(=O)NC2=CC=CC=C2)=C1)C1=CC=CC=C1">(mol);
    validate_contains<"CCCC=C(C1=CC(Cl)=C(OC)C(=C1)C(=O)OC)C1=CC(C(=O)OC)=C(OC)C(Cl)=C1">(mol);
    validate_contains<"CCCCCCCC1CCC(=O)O1">(mol);
    validate_contains<"CCCCCCCCCCC(=C)C(=O)NC1C(OC)C=C(OC)C=C1OC">(mol);
    validate_contains<"CCCCCCCCCCCCCCCC(O)=O">(mol);
    validate_contains<"CCCN1C(N)C(N)C(=O)N(CCC)C1=O">(mol);
    validate_contains<"CCCN1C2NC=NC2C(=O)N(CCC)C1=O">(mol);
    validate_contains<"CCN(C(O)C1CCCC1C(O)NCC1=CC=C(C=C1)C(N)=N)C1CCCCC1">(mol);
    validate_contains<"CCN(CC)CC1=CC=C(CC(NC(=O)NC2=CC(OC)=CC=C2OC)C(O)=O)N1">(mol);
    validate_contains<"CCN(CC)CC1=CC=C(CNC(=O)NC2=CC(OC)=CC=C2OC)O1">(mol);
    validate_contains<"CCN(CC)CCC1=CC=C(CC(NC(=O)NC2=CC(OC)=CC=C2OC)C(O)=O)N1">(mol);
    validate_contains<"CCN(CC)CCC1=CC=C(CNC(=O)NC2=CC(OC)=CC=C2OC)N1">(mol);
    validate_contains<"CCN(CC)CCC1=CC=C(CNC(=O)NC2=CC(OC)=CC=C2OC)O1">(mol);
    //validate_contains<"CCN(CC)CCNC(=O)C1=C(C)NC(C=C2/C(=O)NC3=CC=C(F)C=C23)=C1C">(mol); // FIXME: stereo
    validate_contains<"CCN1C(=O)C(C)=CC2=CN=C(C)N=C12">(mol);
    validate_contains<"CCNC(=O)C1=CC(OC2=CC=C(NC(=O)NC3=CC=CC=C3)C=C2)=CC=N1">(mol);
    validate_contains<"CCNC(=O)C1=NC=CC(OC2=CC=C(NC(N)=O)C=C2)=C1">(mol);
    validate_contains<"CCOC(=O)C1=C(NC(=O)C2=CC=C(Cl)C(Cl)=C2)N=C(SC)N1C1CC=CC=C1">(mol);
    validate_contains<"CCOC(=O)N=S(C)(=O)C1=C2NCCCCNC3=C(Br)C=NC(NC(C=C1)=C2)=N3">(mol);
    //validate_contains<"CCOC(=O)[C@H](CCC1=CC=CC=C1)N[C@@H](C)C(=O)N1[C@H]2CCCC[C@@H]2C[C@H]1C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1Cc2ccccc2C[C@H]1C(O)=O">(mol); // FIXME: stereo
    validate_contains<"CCOc1cc2ncc(C#N)c(Nc3ccc(OCc4ccccn4)c(Cl)c3)c2cc1NC(=O)C=CCN(C)C">(mol);
    validate_contains<"CC[N+](C)(C)C1=CC=CC(O)=C1">(mol);
    validate_contains<"CN(C(=O)CCl)C1=CC=CC=C1">(mol);
    //validate_contains<"CN(C)C(=S)NN=C1/C(=O)NC2=C1C=CC=C2">(mol); // FIXME: stereo
    validate_contains<"CN(C)C1=CC=C(C=C1)C1=[N+](C)C2=C(S1)C=C(O)C=C2">(mol);
    validate_contains<"CN(C)C1=NC=C(C=C1)C1=NC2=NC=NC(N)=C2C(C)=C1">(mol);
    validate_contains<"CN(C)CCC(O)C(N1CCCC1)C1=CC=CC=C1">(mol);
    validate_contains<"CN(C)CCC1=CNC2=C1C(O)=CC=C2">(mol);
}
