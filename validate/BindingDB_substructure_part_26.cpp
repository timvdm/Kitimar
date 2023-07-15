#include "Validate.hpp"

void BindingDB_substructure_part_26(OpenBabel::OBMol &mol)
{
    // SMARTS 1251 - 1300
    validate_contains<"CCC(N)C=CS">(mol);
    validate_contains<"CCC(O)C(O)C(N)CO">(mol);
    //validate_contains<"CCC.CC=C.CC=C.CCC=C/C=CC">(mol); // FIXME: components
    //validate_contains<"CCC.CCCCC1(C)CC1">(mol); // FIXME: components
    //validate_contains<"CCC.OC(=O)C1=CN=C(NC2=CC=NC3=CC(Cl)=CC=C23)C=C1O">(mol); // FIXME: components
    //validate_contains<"CCC.[H]OP(=O)(O[H])O[C@]([H])(C([H])([H])N(C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H])S(=O)(=O)C1=C([H])C([H])=C(N([H])[H])C([H])=C1[H])[C@@]([H])(N([H])C(=O)O[C@]2([H])C([H])([H])OC([H])([H])C2([H])[H])C([H])([H])C3=C([H])C([H])=C([H])C([H])=C3[H]">(mol); // FIXME: components
    validate_contains<"CCC1(CC)CC2=CC=CC=C2C1">(mol);
    validate_contains<"CCC1(CCCC2=C1C=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"CCC1(CN)CCCCC1">(mol);
    validate_contains<"CCC1=C(C)C2=C(OCC(C)=O)C=C(OCC(O)=O)C=C2OC1=O">(mol);
    validate_contains<"CCC1=C(C)NC(=O)C(CC2=NC3=C(O2)C=CC=C3)=C1">(mol);
    validate_contains<"CCC1=C(C)NC(=O)C(N)=C1C(=O)C2=CC(C)=CC(C)=C2">(mol);
    validate_contains<"CCC1=C(C)NC(=O)C(NCN2C(=O)C3=CC=CC=C3C2=O)=C1">(mol);
    validate_contains<"CCC1=C(NN=C1C)C(=O)NC">(mol);
    validate_contains<"CCC1=C2C=CNC2=NC=N1">(mol);
    //validate_contains<"CCC1=C2N=C(C=C(NCC3=CC=C[N+]([O-])=C3)N2N=C1)N4CCCC[C@H]4CCO">(mol); // FIXME: stereo
    validate_contains<"CCC1=CC=C(C(C)C)C2=C1C=C(C=C2)C(C)S(=O)(=O)CC">(mol);
    validate_contains<"CCC1=CC=C(C=C1)S(=O)(=O)NC2=CC(=NN2C3=CC=CC=C3)C4=CC=CS4">(mol);
    validate_contains<"CCC1=CC=C2SC3=C(N=C(CN(C)C)NC3=O)C2=C1">(mol);
    validate_contains<"CCC1=CC=CC2=C1C=CN2">(mol);
    //validate_contains<"CCC1=CC=CC=C1.CCC2=C(C=CC=C2)[N+]([O-])=O">(mol); // FIXME: components
    validate_contains<"CCC1=CC=CN=C1">(mol);
    validate_contains<"CCC1=NC(=C(O)N1C)C2=CC=CC=C2">(mol);
    validate_contains<"CCC1=NNC(C)=N1">(mol);
    validate_contains<"CCC1C2=C(C=CC=C2)C3=C1C=CC=C3">(mol);
    validate_contains<"CCC1CNCCN1">(mol);
    validate_contains<"CCC1OC(O)C(O)C(O)C1OC2OC(CO)C(O)C(O)C2O">(mol);
    validate_contains<"CCC=C">(mol);
    validate_contains<"CCC=NC">(mol);
    validate_contains<"CCCC(C)C1=C2C=C(C)C(C)=CC2=CC3=C1NC4=CC=CC=C34">(mol);
    validate_contains<"CCCC(CC)C(=O)OCCCN">(mol);
    validate_contains<"CCCC(O)=O">(mol);
    validate_contains<"CCCC1=C(C)C2=C(OCC(O)=O)C=C(OCC(O)=O)C=C2OC1=O">(mol);
    validate_contains<"CCCC1=CC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"CCCC1=CC=C(C)C=C1">(mol);
    validate_contains<"CCCC1=CC=C(OCC(=O)C(=CC2=CC=C(OCCO)C=C2)C#N)C=C1">(mol);
    validate_contains<"CCCC1=CC=C2NC3=C(C=C4C=CC=CC4=C3)C2=C1">(mol);
    validate_contains<"CCCC1=CC=C2SC3=C(N=C(CN(C)C)NC3=O)C2=C1">(mol);
    validate_contains<"CCCC1=CN(CC)C2=C1C=CC=C2">(mol);
    validate_contains<"CCCC1=CNC2=C1C=C(Cl)C=C2">(mol);
    validate_contains<"CCCC1=CNC=C1">(mol);
    validate_contains<"CCCC1=NC(CC)=CC(N)=C1C#N">(mol);
    validate_contains<"CCCC1=NCC(=O)N1CC2=CC=C(C=C2)C3=CC=CC=C3">(mol);
    validate_contains<"CCCC1=NN(C)C2=C1N=C(NC2=O)C3=C(OCC)C=CC(=C3)S(=O)(=O)N4CCN(C)CC4">(mol);
    validate_contains<"CCCC1C=CC=C1">(mol);
    validate_contains<"CCCC1CCCN(C)C1">(mol);
    validate_contains<"CCCC2(CCC1=CC=CC=C1)CC(=O)C(C(CC)C3=CC=CC(NS(=O)(=O)C4=NC=C(C=C4)C(F)(F)F)=C3)=C(O)O2">(mol);
    validate_contains<"CCCC=C(C1=CC(C)=C(OC)C(=C1)C(=O)OC)C1=CC(C(=O)OC)=C(OC)C(C)=C1">(mol);
    validate_contains<"CCCCC(=O)NO">(mol);
    validate_contains<"CCCCC(C)CCC">(mol);
}
