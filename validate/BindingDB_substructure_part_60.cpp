#include "Validate.hpp"

void BindingDB_substructure_part_60(OpenBabel::OBMol &mol)
{
    // SMARTS 2951 - 3000
    validate_contains<"CNCC1=CNC=C1">(mol);
    validate_contains<"CNCCCC">(mol);
    validate_contains<"CNS(=O)(=O)c1ccc(Nc2nccc(n2)-c3cnc4ccccn34)cc1">(mol);
    //validate_contains<"CN[C@H](C(c1ccccc1)c2ccccc2)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CCCNC(N)=N)C(=O)c4nc5ccccc5s4">(mol); // FIXME: stereo
    validate_contains<"COC">(mol);
    validate_contains<"COC(=O)CNC(=O)C1=CC2=CC=CC=C2OC1=O">(mol);
    //validate_contains<"COC1=C(O)C=CC(C=C2/CCCC2=O)=C1">(mol); // FIXME: stereo
    validate_contains<"COC1=C(OC)C=C2CN(C)CCC2=C1">(mol);
    validate_contains<"COC1=C(OCCCN2CCCCC2)C=C3N=CN=C(NC4=CN=C(NC(=O)C5=CC=CC=C5)N=C4)C3=C1">(mol);
    validate_contains<"COC1=CC(C)=CC(OC)=C1OC">(mol);
    validate_contains<"COC1=CC(NC2=C(C=NC3=CC(OC)=C(OC)C=C23)C#N)=C(Cl)C=C1Cl">(mol);
    validate_contains<"COC1=CC2=C(C=C1)C3=C(N2)C(C)=NC=C3">(mol);
    validate_contains<"COC1=CC=NC=C1OC">(mol);
    validate_contains<"COCC1(C)CC2=C(C)C(N)=C(C)C(C)=C2O1">(mol);
    //validate_contains<"COC[C@@H]1CCCN1S(=O)(=O)c2ccc3NC(=O)C(=O)c3c2">(mol); // FIXME: stereo
    validate_contains<"COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O">(mol);
    validate_contains<"COc1cc2CC(CC3CCN(CC3)Cc4ccccc4)C(=O)c2cc1OC">(mol);
    validate_contains<"COc1cc2N(CCO)C=C(C(O)=O)C(=O)c2cc1Cc3cccc(Cl)c3F">(mol);
    validate_contains<"COc1cc2c(Nc3ccc(Cl)cc3F)ncnc2cc1OCCn4cncn4">(mol);
    validate_contains<"COc1cc2c(Oc3ccc(nc3)C4=CN=C(Nc5ccc(F)cc5)N(C)C4=O)ccnc2cc1OCCCN6CCOCC6">(mol);
    //validate_contains<"COc1ccc(cc1)-c2[nH]c(nc2CCNS(=O)(=O)N[C@@H]3CCN(C3)Cc4ccccc4)-c5ccccc5">(mol); // FIXME: stereo
    //validate_contains<"COc1ccc2CN(C)CC[C@@]34C=C[C@H](O)C[C@@H]3Oc1c24">(mol); // FIXME: stereo
    validate_contains<"COc1ccc2[nH]c(I)c(CCNC(C)=O)c2c1">(mol);
    validate_contains<"COc1ccc2cc(oc2c1)C(c3ccc(C)cc3)n4cncn4">(mol);
    validate_contains<"CSC1=NC(C2=CC=CC=C2)=C(N=N1)C3=CC=CC=C3">(mol);
    validate_contains<"CSC1=NC=CC=N1">(mol);
    //validate_contains<"CSCC[C@H](N)C(=O)NC1=CC=C2C(C)=CC(=O)OC2=C1">(mol); // FIXME: stereo
    //validate_contains<"CSC[C@H](NC(=O)COc1cccc2cnccc12)C(=O)N[C@@H](Cc3ccccc3)[C@H](O)C(=O)N4CSC[C@H]4C(=O)NC(C)(C)C">(mol); // FIXME: stereo
    //validate_contains<"C[C@@]12CCC[C@@H]1C3CCC4CCCCC4C3CC2">(mol); // FIXME: stereo
    //validate_contains<"C[C@H](NC(=O)[C@H](Cc1ccc(OP(O)(O)=O)cc1)NC(C)=O)c2nc(Cc3ccc(F)cc3)no2">(mol); // FIXME: stereo
    //validate_contains<"C[C@H](NC1=NC(=O)[C@](C)(S1)C(F)(F)F)c2ccccc2F">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1CCCC[C@@H]1Oc2cccc(c2O)-c3cc4cc(C(N)=N)c(F)cc4[nH]3">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1Cn2c(S)nc3cc(Cl)cc(CN1CC=C(C)C)c23">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n2cc(I)c3c(N)ncnc23">(mol); // FIXME: stereo
    //validate_contains<"C[C@](O)([C@H]1CCCC2=Cc3c(C[C@]12C)cnn3-c4ccc(F)cc4)c5ccsc5">(mol); // FIXME: stereo
    //validate_contains<"C[C@]12CCC[C@H]1C3CCC4CCCCC4C3CC2">(mol); // FIXME: stereo
    //validate_contains<"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC(CCCCCCCC(=O)OCC5OC(C(O)C5O)n6cnc7c(N)ncnc67)[C@@H]2O">(mol); // FIXME: stereo
    //validate_contains<"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CCC2=O">(mol); // FIXME: stereo
    validate_contains<"Cc1ccc(F)c(NC(=O)Nc2ccc(cc2)-c3cccc4[nH]nc(N)c34)c1">(mol);
    validate_contains<"Cc1cccc2n(Cc3c(F)cccc3F)c(nc12)-c4c(F)cccc4F">(mol);
    validate_contains<"Cc1nccs1">(mol);
    //validate_contains<"Cl.CC1(C)CCN(CC1)c2ccc(O)c(c2)C(=O)c3ccc(cc3)C(=O)N[C@@H]4CCCNC[C@H]4NC(=O)c5ccncc5">(mol); // FIXME: components
    //validate_contains<"Cl.CCOC1=C(OCC2=C(Cl)C=CC=C2)C=C(Br)C(CN)=C1">(mol); // FIXME: components
    //validate_contains<"Cl.N[C@@H](CC(=O)N1CCn2cnnc2C1)Cc3cc(F)c(F)cc3F">(mol); // FIXME: components
    validate_contains<"ClC1=C(Cl)C=C(COC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"ClC1=CC=C(NC=O)N=C1">(mol);
    validate_contains<"ClC1=NC=CS1">(mol);
    validate_contains<"ClC1=NC=NC=C1">(mol);
    validate_contains<"ClC1=NNC(=C1)C1=NC2=C(C=CC=C2)N2C(=NN=C12)C1=CC=CC=C1">(mol);
    validate_contains<"Clc1ccccc1Nc2ncnc3ccccc23">(mol);
}
