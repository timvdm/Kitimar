#include "Validate.hpp"

void BindingDB_substructure_part_20(OpenBabel::OBMol &mol)
{
    // SMARTS 951 - 1000
    validate_contains<"CC(N)C(C)=O">(mol);
    validate_contains<"CC(N)C(C=C)C1=CC=CC=C1">(mol);
    validate_contains<"CC(N)C(O)=O">(mol);
    validate_contains<"CC(N)C(O)NCC(O)=O">(mol);
    validate_contains<"CC(N)C([O-])=O">(mol);
    validate_contains<"CC(N)CC1=CC=CC(F)=C1">(mol);
    validate_contains<"CC(N)CNC1=CC=CC=C1">(mol);
    validate_contains<"CC(NC(=O)C1CC2=CC=CC=C2CN1C)C(O)=O">(mol);
    validate_contains<"CC(NC(=O)C1CCCCN1C)C(O)=O">(mol);
    validate_contains<"CC(NC(=O)NC(Cl)(Cl)Cl)C1=CC(=O)C2=CC=CC=C2O1">(mol);
    validate_contains<"CC(NC1=NC(=O)C(S1)C(F)(F)F)C2=C(F)C=CC=C2">(mol);
    validate_contains<"CC(NS)C=O">(mol);
    validate_contains<"CC(O)(CSC1=CC=C(C=C1)[N+](O)=O)C(=O)NC2=CC(=C(C=C2)C#N)C(F)(F)F">(mol);
    validate_contains<"CC(O)C(C)=O">(mol);
    validate_contains<"CC(O)C(C)C(N)C(O)=O">(mol);
    validate_contains<"CC(O)C(NC(=O)C1CCCN1C(=O)C(N)CCCCN)C(O)=O">(mol);
    validate_contains<"CC(O)C1=CC(O)=CC(O)=C1">(mol);
    validate_contains<"CC(O)C1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CC(O)CC(C)(C)C(CC#C)C(C)(C)C">(mol);
    validate_contains<"CC(O)CC(C)(C)C(CCC#C)C(C)(C)C">(mol);
    validate_contains<"CC(O)CC1=CC2=C(OCO2)C=C1[N+](O)=O">(mol);
    validate_contains<"CC(O)NC1=CC=CC=C1">(mol);
    validate_contains<"CC(O)OC(=O)CCCCC1CCSS1">(mol);
    validate_contains<"CC(OC1=C(N)N=CC(=C1)C2=CC=CC=C2)C3=C(Cl)C(F)=CC=C3Cl">(mol);
    //validate_contains<"CC.C(C1CCCCC1)N2CC=CC2">(mol); // FIXME: components
    //validate_contains<"CC.CC(N)CN(C[C@@H](O)[C@H](CC1=CC(O)=CC=C1)NC(O)=O)S(=O)(=O)C2=CC=C(O)C=C2">(mol); // FIXME: components
    //validate_contains<"CC.CC.CC">(mol); // FIXME: components
    //validate_contains<"CC.CC.CC.C1CCCC1.C2CCCC2.C3CCCC3.C4CCCC4.C5CCCC5">(mol); // FIXME: components
    //validate_contains<"CC.CC.CC.CC.C1CCCC1">(mol); // FIXME: components
    //validate_contains<"CC.CC.CC.CC.CC.CC.CCC.CCC">(mol); // FIXME: components
    //validate_contains<"CC.CC.CC.CC.CC.CCC">(mol); // FIXME: components
    //validate_contains<"CC.CC.CC1=CC=C(O)C=C1.CC2=CC(O)=CC(O)=C2">(mol); // FIXME: components
    //validate_contains<"CC.CC.CCC">(mol); // FIXME: components
    //validate_contains<"CC.CC1=C(N2C(SC1)C(NC(=O)CCCC(N)C(O)=O)C2=O)C(O)=O">(mol); // FIXME: components
    //validate_contains<"CC.CC1=CC=CC=C1">(mol); // FIXME: components
    //validate_contains<"CC.CC1=CC=CC=C1CNC(=O)[C@H]2N(CSC2(C)C)C(=O)[C@@H](O)[C@H](CC3=CC=CC=C3)NC(=O)C4=C(C)C(O)=CC=C4">(mol); // FIXME: components
    //validate_contains<"CC.CC1CCC(C)C1.CC2CCC(C)C2">(mol); // FIXME: components
    //validate_contains<"CC.CCC(CC)CC">(mol); // FIXME: components
    //validate_contains<"CC.CCOC1OC(C)(C2CCC3C4CC5OC56C(O)C=CC(=O)C6C4CCC123)C7CC(C)=C(C)C(=O)O7">(mol); // FIXME: components
    //validate_contains<"CC.CNC.N1C=CC=C1.N2C=CC=C2.C3CC4CCC5(CCCC5)C4C3.C6CC7=CC=CC8=C7C6=CC=C8">(mol); // FIXME: components
    //validate_contains<"CC.NC(N)=N/C1=SC=CN1">(mol); // FIXME: components
    //validate_contains<"CC.O[C@@H]1CCC(N2C(=O)NC3=CN=C(N=C23)N4C=NC5=C4C=C(C=C5)C#N)C6=C1C(F)=CC(F)=C6">(mol); // FIXME: components
    //validate_contains<"CC.O[C@@H]1[C@@H](O)[C@@H](CC2=CC=CC=C2)N3CC4=CC=C5C=CC=CC5=C4C6=CC(CN([C@@H]1CC7=CC=CC=C7)C3=O)=CC8=CC=CC=C68">(mol); // FIXME: components
    //validate_contains<"CC.[H]N(N([H])C1=NC(=NC(=C1[H])C([H])([H])[H])N([H])S(=O)(=O)C2=C([H])C([H])=C([H])C([H])=C2[H])C([H])=C3/C([H])=C([H])C(=O)C(OC([H])([H])[H])=C3[H]">(mol); // FIXME: components
    //validate_contains<"CC.[H]N1C(C)C(C)(C)N([H])C1=O">(mol); // FIXME: components
    //validate_contains<"CC.c1ccc2c(c1)ccc3c4ccccc4ccc23">(mol); // FIXME: components
    validate_contains<"CC1(C)C(=O)NC(SC2CCCCC2)=NC1=CC3=CC=CC=C3">(mol);
    validate_contains<"CC1(C)C(C(=O)C2=CN(CCN3CCOCC3)C4=CC=CC=C24)C1(C)C">(mol);
    validate_contains<"CC1(C)CC2=CC=CC=C2C1">(mol);
    validate_contains<"CC1(C)CCC(C)(C)C2=C1C=CC=C2">(mol);
}
