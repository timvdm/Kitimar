#include "Validate.hpp"

void BindingDB_substructure_part_58(OpenBabel::OBMol &mol)
{
    // SMARTS 2851 - 2900
    //validate_contains<"CC(C)C[C@H](NC(=O)[C@@H](N)C(C)C)C(=O)N1CCCC1C(=O)N[C@@H](CCC2=CC=CC=C2)C(N)=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)NCC(O)c1ccc(O)c(O)c1">(mol);
    //validate_contains<"CC(C)NC[C@H](O)COc1ccc(CCOCC2CC2)cc1">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)N(C)Cc1csc(n1)C(C)C)C(=O)N[C@H](C[C@H](O)[C@H](Cc2ccccc2)NC(=O)OCc3cncs3)Cc4ccccc4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(C)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(O)=O)C(=O)CCCCCc2ccccc2">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(C)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(O)=O)C=O">(mol); // FIXME: stereo
    //validate_contains<"CC(C)c1c(C(=O)Nc2ccccc2)c(-c3ccccc3)c(-c4ccc(F)cc4)n1CC[C@@H](O)C[C@@H](O)CC(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(C)c1nc(nc(-c2ccc(F)cc2)c1C=C[C@@H](O)C[C@@H](O)CC(O)=O)N(C)S(C)(=O)=O">(mol); // FIXME: stereo
    validate_contains<"CC(C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CC(N)CN1C=CC=CC1=O">(mol);
    validate_contains<"CC(O)C(C)NC(=O)c1cc(cc(c1)C(=O)NC(C)c2ccccc2)N(C)S(C)(=O)=O">(mol);
    //validate_contains<"CC.CC.CC.CC.CC.CC.CC.CC.CC.CC.CC.CCC.CCC.CC=C.CCCC.CCCCCC.CCCCCC">(mol); // FIXME: components
    validate_contains<"CC1(C)CC(=S)NC2=CC=C(C=C12)C3=CC=C(N3)C#N">(mol);
    validate_contains<"CC1(CCCNCC1)c2nc3c(cccc3[nH]2)C(N)=O">(mol);
    //validate_contains<"CC1(CCN(CC1)c2cc(ccn2)C(O)=O)NCC(=O)N3[C@H](CC[C@H]3C#N)C#C">(mol); // FIXME: stereo
    validate_contains<"CC1=C(SC(N)=N1)C(=O)NC2=CC=CC=C2">(mol);
    validate_contains<"CC1=CC(NC2=NC=CC(=N2)N3C=CN=C3C4=CC=CC=C4)=CC(C)=C1">(mol);
    validate_contains<"CC1=CC(NS(=O)(=O)C2=CC=CC=C2)=CC(OCC[NH2+]C3=C=CNC=C3)=C1">(mol);
    validate_contains<"CC1=CC2=C(C=C1)C(Cl)=CC(F)=C2">(mol);
    validate_contains<"CC1=CC2=C(NC=C2)C=C1">(mol);
    validate_contains<"CC1=CC2=C(OC(=O)C(=C2)C(=O)NC3=CC(C)=CC(C)=C3)C=C1">(mol);
    //validate_contains<"CC1=CC=CC(CNC(=O)[C@@H]2CCCN2C(=O)[C@H]([NH3+])CC3=CC=CC=C3)=C1">(mol); // FIXME: stereo
    validate_contains<"CC1=CC=CC2=C1C=CC=C2">(mol);
    //validate_contains<"CC1=CN([C@H]2C[C@H](O)[C@@H](CO)O2)C(=O)NC1=O">(mol); // FIXME: stereo
    validate_contains<"CC1=CSC(N)=N1">(mol);
    validate_contains<"CC1=NC(C)=NC(C)=N1">(mol);
    //validate_contains<"CC1=NC(CC(=N/O)C2=C(O)C(O)=C(O)C=C2)=CS1">(mol); // FIXME: stereo
    validate_contains<"CC1=NC(N)=NC=C1">(mol);
    validate_contains<"CC1=NOC(=N1)C1=CNC2=CC=CC=C12">(mol);
    validate_contains<"CC1C2=CC=CC=C2CCC3=C1N=CC=C3">(mol);
    validate_contains<"CC1CCN(C)CC1O">(mol);
    validate_contains<"CC1CN(C)C(=O)O1">(mol);
    validate_contains<"CC1CNCC1C">(mol);
    validate_contains<"CC=CC1=CC=C(C=C1)S(N)(=O)=O">(mol);
    //validate_contains<"CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)OCCCc2cccnc2">(mol); // FIXME: stereo
    validate_contains<"CCC(C)C">(mol);
    validate_contains<"CCC(C1=CC=C(O)C=C1)=C(CC)C2=CC=C(O)C=C2">(mol);
    validate_contains<"CCC(Cc1ccccc1)NC=O">(mol);
    validate_contains<"CCCC(=O)NO">(mol);
    validate_contains<"CCCC1=CC=C(CC(=C)C2=CC=CC=C2CCC3=CNC(=O)C=C3)C=C1C(F)(F)F">(mol);
    validate_contains<"CCCC=C1SC2=C(NC3=CC=CC=C23)cc1">(mol);
    validate_contains<"CCCCCCCCOC(=O)C1=CC=CC=C1C(=O)OCCCC">(mol);
    validate_contains<"CCCCCCO">(mol);
    validate_contains<"CCCCCCc1cc2C(O)C(C(c3ccccc3)C4=C(O)c5cc(CCCCCC)c(O)cc5OC4=O)C(=O)Oc2cc1O">(mol);
    validate_contains<"CCCl">(mol);
    validate_contains<"CCN1C=CC2=C1CCCC2=O">(mol);
    validate_contains<"CCN1c2nc(ccc2N(C)C(=O)c3cccnc13)-c4cccc(N)c4">(mol);
    validate_contains<"CCOC(=O)C1=CN(CC2=CC=CC=C2)C=C(C1C3=CC=CC=C3)C(=O)OCC">(mol);
    validate_contains<"CCOC(=O)C1=NOC=C1">(mol);
    validate_contains<"CCOC(=O)C=CC(CC1CCNC1=O)NC(=O)C(CC(=O)C(NC(=O)C2=NOC(C)=C2)C(C)C)CC3=CC=C(F)C=C3">(mol);
}
