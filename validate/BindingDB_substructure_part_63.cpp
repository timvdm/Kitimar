#include "Validate.hpp"

void BindingDB_substructure_part_63(OpenBabel::OBMol &mol)
{
    // SMARTS 3101 - 3150
    validate_contains<"OC(=O)C1=CC(=CC=C1)C(O)=O">(mol);
    validate_contains<"OC(=O)C1=CC=C(NC=O)C=C1">(mol);
    validate_contains<"OC(=O)C1CC1c2ccc(NCc3ccc(cc3)-c4ccccc4)cc2">(mol);
    validate_contains<"OC(=O)Cc1c(SSc2[nH]c3ccccc3c2CC(O)=O)[nH]c4ccccc14">(mol);
    validate_contains<"OC(O)C1CCCO1">(mol);
    validate_contains<"OC1=C(CCCOC2CCCCC2)C(=O)OC3=C1C=CC=C3">(mol);
    validate_contains<"OC1=C(N=C2C(CCCN2C1=O)N3CCOCC3)C(=O)NCc4ccc(F)cc4">(mol);
    validate_contains<"OC1=CC2=CC=CC=C2C(O)=C1">(mol);
    validate_contains<"OC1=CC=C2C(OC(=O)C=C2C(=O)NC3=CC=CC=C3)=C1">(mol);
    validate_contains<"OC1=CC=CC2=C1C(O)=C(CC3=CC=CC=C3)C(=O)O2">(mol);
    validate_contains<"OC1=CC=CC2=C1NCCC2">(mol);
    validate_contains<"OC1=CC=NC2=C1C=CC(Cl)=C2">(mol);
    validate_contains<"OC1=CC=NC2=CC=CC=C12">(mol);
    validate_contains<"OC1CCC(=O)C1">(mol);
    validate_contains<"OC1CCC(O)C1">(mol);
    validate_contains<"OC1CCC2C1CCC3C2CCC4=C3C=CC(O)=C4">(mol);
    validate_contains<"OC1CCC2C1CCC3C2CCC4=CC(O)=CC=C34">(mol);
    validate_contains<"OC1CNC=NC2=C1N=CN2C3CCCO3">(mol);
    validate_contains<"OCC(NC(=O)C1CCCN1C=O)C(O)=O">(mol);
    validate_contains<"OCC(NC(=O)c1cc(c[nH]1)-c2n[nH]cc2-c3cccc(Cl)c3)c4ccc(Cl)cc4">(mol);
    validate_contains<"OCC1CC(O)C(O)CO1">(mol);
    validate_contains<"OCCCC#CCCCCCC(O)=O">(mol);
    validate_contains<"ON=C1C(Nc2cc(Br)ccc12)=C3C(=O)Nc4cc(Br)ccc34">(mol);
    validate_contains<"Oc1cc(O)c2C(=O)C(O)=C(Oc2c1)c3cc(O)c(O)c(O)c3">(mol);
    validate_contains<"Oc1cc2C(=O)Oc3c(O)c(O)cc4C(=O)Oc(c1O)c2-c34">(mol);
    validate_contains<"Oc1ccc(cc1)-c2ccc3cc(O)ccc3c2Cc4ccc(OCCN5CCCCC5)cc4">(mol);
    validate_contains<"Oc1ccc2[nH]c(cc2c1)C(=O)c3cc4ccccc4[nH]3">(mol);
    validate_contains<"S1C=NC(=C1)N1C=CC=N1">(mol);
    validate_contains<"SC1=CNC2=C1N=CC=C2">(mol);
    validate_contains<"SC1=NC=CC2=C1C=CC=C2">(mol);
    validate_contains<"SC1=NC=CC=C1Br">(mol);
    validate_contains<"SCCCC(=O)NC1CCN(CC2=CC=CC=C2)CC1">(mol);
    //validate_contains<"[Br-].C[n+]1c(-c2ccccc2)c3cc(N)ccc3c4ccc(N)cc14">(mol); // FIXME: components
    //validate_contains<"[Fe].[cH-]1[cH-][cH-][cH-][cH-]1.OC(=O)C[c-]2cccc2">(mol); // FIXME: components
    validate_contains<"[H]C(=O)C(O)COP([O-])([O-])=O">(mol);
    validate_contains<"[H]C12CCCCC1([H])CN(CC(O)C(CSC3=CC=CC=C3)NC(=O)C4=CC=CC(O)=C4C)C(C2)C(=O)NC(C)(C)C">(mol);
    validate_contains<"[H]C1=C([H])C(=NC(NC)=N1)C2=CC=CC=C2">(mol);
    validate_contains<"[H]C1=NN(*)C([H])=N1">(mol);
    validate_contains<"[H]N(C1=NC(C)=C([H])S1)C2=CC=CC=C2">(mol);
    validate_contains<"[H]N([H])C(NC)=NC">(mol);
    validate_contains<"[H]N([H])C1=NC2=C(N=C([H])N2[H])C(OC([H])([H])C3([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C3([H])[H])=N1">(mol);
    validate_contains<"[H]N([H])S(=O)(=O)C1=CC=CC=C1">(mol);
    validate_contains<"[H]N1CC2=C3C4=CC=CC=C4[N]C3=C5[N]C6=CC=CC=C6C5=C2C1=O">(mol);
    validate_contains<"[H]OC">(mol);
    validate_contains<"[H]OC1=C([H])C(SC(=C([H])C(=O)OC([H])([H])C([H])([H])[H])C([H])([H])[H])=C(O[H])C2=C1C([H])=C([H])C([H])=C2[H]">(mol);
    validate_contains<"[H]OCC[N+]([H])([H])C">(mol);
    //validate_contains<"[H][C@@](CCCNC(N)=[NH2+])(NS(=O)(=O)C1=CC=CC2=C1C=CC=C2N(C)C)C(=O)N3CCCC[C@]3([H])CCC(=O)C4=NC=CS4">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CCCC[C@]1([H])CN(C[C@@H](O)[C@H](CC3=CC=CC=C3)NC(=O)[C@H](CC(N)=O)NC(=O)C4=NC5=CC=CC=C5C=C4)[C@@H](C2)C(=O)NC(C)(C)C">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CCC[C@@]1(C)CCC3C4CCCCC4CCC23">(mol); // FIXME: stereo
    //validate_contains<"[I-].CC(=O)SCC[N+](C)(C)C">(mol); // FIXME: components
}
