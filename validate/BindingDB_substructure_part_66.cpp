#include "Validate.hpp"

void BindingDB_substructure_part_66(OpenBabel::OBMol &mol)
{
    // SMARTS 3251 - 3300
    validate_contains<"c1ccc2ccccc2c1">(mol);
    validate_contains<"CC1=NCCCN1">(mol);
    validate_contains<"N1C=CN=C1">(mol);
    validate_contains<"NC1=CSC=N1">(mol);
    validate_contains<"C1=CC2=C(C=C1)N=CN=C2">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"C1CC2=CC=CC=C2N1">(mol);
    validate_contains<"C1CC=NN1">(mol);
    validate_contains<"C1COC2=CC=CC=C2O1">(mol);
    validate_contains<"CC(C)CN(CC(O)C(CC1=CC=CC=C1)NC(=O)OC2CCOC2)S(=O)(=O)C3=CC=C(N)C=C3">(mol);
    //validate_contains<"CC.CC">(mol); // FIXME: components
    validate_contains<"CC12CCCC1C3CCC4=CCCCC4C3CC2">(mol);
    validate_contains<"CC1=CC=C(Cl)S1">(mol);
    validate_contains<"CC1=CC=CN=C1N">(mol);
    validate_contains<"CC1CCCN1">(mol);
    validate_contains<"CCC(CO)NC1=NC2=C(N=CN2C(C)C)C(NCC3=CC=CC=C3)=N1">(mol);
    validate_contains<"CCCCCC=CC=CC=O">(mol);
    validate_contains<"CCCCCCCCCCCCCCCCCCCCCCN1CCN(CC1)C(=O)c2ccc(CC3=NOC(=O)N3)cc2">(mol);
    //validate_contains<"CC[C@@H]">(mol); // FIXME: stereo
    validate_contains<"CN(C)c1ccc(cc1)C#Cc2ncnc(N)c2-c3ccc4OCOc4c3">(mol);
    validate_contains<"CN1CC=CN=C1">(mol);
    validate_contains<"CN1CCCC1C2=CC=CN=C2">(mol);
    validate_contains<"CN1CCN(CC1)c2cc(Nc3cc(C)n[nH]3)nc(Sc4ccc(NC(=O)C5CC5)cc4)n2">(mol);
    //validate_contains<"CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)c3ccc(C)cc3)c4ncc(s4)-c5cccc(c5)[N+]([O-])=O">(mol); // FIXME: stereo
    validate_contains<"CNC(C)=O">(mol);
    validate_contains<"COC1=C(O)C=C2N=CN=C(N)C2=C1">(mol);
    validate_contains<"COP">(mol);
    validate_contains<"COc1cc(NC(=S)NCCCn2ccnc2)cc(OC)c1OC">(mol);
    //validate_contains<"C[C@@H](O)[C@H](N)CC1=CC=CC=C1">(mol); // FIXME: stereo
    //validate_contains<"Cc1c(O)cccc1C(=O)N[C@@H](CSc2ccccc2)[C@H](O)CN3C[C@H]4CCCC[C@H]4C[C@H]3C(=O)NC(C)(C)C">(mol); // FIXME: stereo
    //validate_contains<"Clc1cccc(Cl)c1N=C2/NCCN2">(mol); // FIXME: stereo
    validate_contains<"N(C1=CC=CC=C1)C2=NC=CC=N2">(mol);
    validate_contains<"N1C2=CC=CC=C2C3=C1C=CC=C3">(mol);
    //validate_contains<"N1C=CC=C1.C2=CC=C3C=CC=CC3=C2">(mol); // FIXME: components
    //validate_contains<"N1C=CC=C1.C2=CC=CC=C2.C3=CC4=C(C=C3)C=CC=C4">(mol); // FIXME: components
    validate_contains<"N1N=CC2=CN=CN=C12">(mol);
    validate_contains<"NC(=N)C1=CC=C(C=C)C=C1">(mol);
    validate_contains<"NC(=N)c1ccc2nc(Cc3nc4ccc(cc4[nH]3)C(N)=N)[nH]c2c1">(mol);
    validate_contains<"NC(=O)C1=C(N)C=CC=C1">(mol);
    //validate_contains<"NC(=[NH2+])C1=CC=C(CNC(=O)[C@@H]2CCCN2C(=O)CCC3=CC=CC=C3)C=C1">(mol); // FIXME: stereo
    validate_contains<"NC1=C2CCCCC2=NC3=C1C=CC=C3">(mol);
    validate_contains<"NC1=CC=C(C=C1)S(N)(=O)=O">(mol);
    validate_contains<"NC1=NC(N)=NC=C1">(mol);
    validate_contains<"NC1=NC=CC=C1">(mol);
    validate_contains<"NC1=NC=NN1">(mol);
    validate_contains<"NC1=[S]C2=C(N1)C(Cl)=CC=C2">(mol);
    validate_contains<"NCCC1=CNC=N1">(mol);
    validate_contains<"O(C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CC2=CN=CN=C2N1">(mol);
    validate_contains<"O=C1CCCO1">(mol);
}
