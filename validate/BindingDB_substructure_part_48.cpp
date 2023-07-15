#include "Validate.hpp"

void BindingDB_substructure_part_48(OpenBabel::OBMol &mol)
{
    // SMARTS 2351 - 2400
    validate_contains<"O=CNC1=NCCS1">(mol);
    validate_contains<"O=N(=O)C1CCCC1">(mol);
    validate_contains<"O=S(=O)(N1CCCC1)c2ccc(CCc3ccc(OCc4ccccc4)cc3)cc2">(mol);
    validate_contains<"O=S(=O)(N1CCCC2=CC=CC=C12)C3=CC=CC=C3">(mol);
    validate_contains<"O=S(=O)(N1CCCNCC1)C2=CC=CC3=CN=CC=C23">(mol);
    validate_contains<"O=S(=O)(NC1=CC=CC=C1)C1=C2C=CN=CC2=CC=C1">(mol);
    validate_contains<"O=S(=O)(NC1=NC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"O=S(=O)(NC1=NC=CS1)C2=CC=CC=C2">(mol);
    validate_contains<"O=S(=O)(NCC1CCNCC1)C2=CC=C(S2)C3=CC=CC=C3">(mol);
    validate_contains<"O=S1(=O)N=CNC2=C1C=CC=C2">(mol);
    validate_contains<"O=S1(=O)NC=NC2=C1C=CC=C2">(mol);
    validate_contains<"O=S1(=O)NCNC2=CC=CC=C12">(mol);
    //validate_contains<"OC(=C/C(=O)c1ccc(Cc2ccc(F)cc2)o1)c3nc[nH]n3">(mol); // FIXME: stereo
    validate_contains<"OC(=O)C(=O)CC(=O)c1cccn1Cc2ccc(F)cc2">(mol);
    //validate_contains<"OC(=O)C(=O)CC1=CC=C(O)C=C1.OC(=O)C(=O)CC2=CNC3=C2C=CC=C3">(mol); // FIXME: components
    validate_contains<"OC(=O)C(=O)CC1=CNC2=C1C=CC=C2">(mol);
    validate_contains<"OC(=O)C(=O)NC1=CC=CC=C1C(O)=O">(mol);
    validate_contains<"OC(=O)C(C(O)=O)C1=CC=CC=C1">(mol);
    validate_contains<"OC(=O)C(CC1=CC=CC=C1)OC2=C(Br)C=C(C=C2Br)C3=C4C5=C(OC4=C(Br)C6=CC=CC=C36)C=CC=C5">(mol);
    validate_contains<"OC(=O)C(O)=CC=O">(mol);
    //validate_contains<"OC(=O)C1=C(C(O)=CC=C1)C(=O)C2=C(O)C=C(C=C2O)C(=O)O[C@@H]3C4CC(C5CCCC45)[C@H]3NC(=O)C6=CC=C(O)C=C6">(mol); // FIXME: stereo
    //validate_contains<"OC(=O)C1=C(C(O)=CC=C1)C(=O)C2=C(O)C=C(C=C2O)C(=O)O[C@@H]3[C@@H]4C[C@@H]([C@H]5CCC[C@@H]45)[C@H]3NC(=O)C6=CC=C(O)C=C6">(mol); // FIXME: stereo
    validate_contains<"OC(=O)C1=C(C=CC(=C1)C(=O)NCC2=CC=C(C=C2)C3=CNN=C3C4=CC(Cl)=C(O)C=C4O)C5=C6C=CC(=O)C=C6OC7=C5C=CC(O)=C7">(mol);
    validate_contains<"OC(=O)C1=C(C=CC=C1)C(O)=O">(mol);
    validate_contains<"OC(=O)C1=C(C=CC=C1)[N+]([O-])=O">(mol);
    validate_contains<"OC(=O)C1=C(O)C2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"OC(=O)C1=C2N=CC=NC2=CC=C1">(mol);
    validate_contains<"OC(=O)C1=CC(=O)OC2=C1C=CC(O)=C2">(mol);
    validate_contains<"OC(=O)C1=CC(=O)OC2=C1C=CC=C2">(mol);
    validate_contains<"OC(=O)C1=CC(CN2CCCC2)=CC=C1">(mol);
    validate_contains<"OC(=O)C1=CC2=CC3=C(C=CC=C3)C=C2OC1=O">(mol);
    validate_contains<"OC(=O)C1=CC2=CC=CC=C2C=C1">(mol);
    validate_contains<"OC(=O)C1=CC=C(C=C1)C2=CC=C3C=CC=CC3=C2">(mol);
    //validate_contains<"OC(=O)C1=CC=CC(O)=C1C(=O)C2=C(O)C=C(C=C2O)C(=O)O[C@@H]3[C@@H]4C[C@@H](C5CCCC45)[C@H]3NC(=O)C6=CC=C(O)C=C6">(mol); // FIXME: stereo
    //validate_contains<"OC(=O)C1=CC=CC(O)=C1C(=O)C2=C(O)C=C(C=C2O)C(=O)O[C@@H]3[C@H]4C[C@H]([C@H]5CCC[C@@H]45)[C@H]3NC(=O)C6=CC=C(O)C=C6">(mol); // FIXME: stereo
    validate_contains<"OC(=O)C1=CC=CC2=C1N=C(O2)C3=CC=CC=C3">(mol);
    validate_contains<"OC(=O)C1=CC=CC=C1CNC=O">(mol);
    validate_contains<"OC(=O)C1=CNc2ccc(Cc3cccc(Cl)c3F)cc2C1=O">(mol);
    validate_contains<"OC(=O)C1=NOC(=C1)C2=CC=CC=C2O">(mol);
    validate_contains<"OC(=O)C1CC2C(NC3=CC=CC=C23)C=N1">(mol);
    validate_contains<"OC(=O)C1CCCC1">(mol);
    validate_contains<"OC(=O)C1CNC(C1)C(O)=O">(mol);
    validate_contains<"OC(=O)C1CNC2CC3=CNC4=CC=CC(C2C1)=C34">(mol);
    validate_contains<"OC(=O)C1CNCCN1">(mol);
    validate_contains<"OC(=O)C=Cc1ccc(O)c(OC2=Cc3cc(O)c(O)cc3OC2=O)c1">(mol);
    validate_contains<"OC(=O)C=Cc1ccc(OC2=Cc3cc(O)c(O)cc3OC2=O)c(O)c1">(mol);
    validate_contains<"OC(=O)CC1=CC(C2=CNC3=C2C=CC=C3)=C(O)C(=C1)C4=CC=CC=C4">(mol);
    validate_contains<"OC(=O)CC1=CC=C(COC2=CC=CC(=C2)C3=C(CC4=CC=CC=C4)C=NC5=C3C=CC=C5C(F)(F)F)C=C1">(mol);
    validate_contains<"OC(=O)CCC(=O)CSCC1=CC=CC=C1">(mol);
    validate_contains<"OC(=O)CCC(=O)NC1=CC=C(C=C1)S(=O)(=O)C2=CC=C(NC(=O)CCC(O)=O)C=C2">(mol);
}
