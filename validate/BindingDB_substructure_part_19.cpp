#include "Validate.hpp"

void BindingDB_substructure_part_19(OpenBabel::OBMol &mol)
{
    // SMARTS 901 - 950
    validate_contains<"CC(C)C1=CC=CC(=C1)C(N)=N">(mol);
    validate_contains<"CC(C)CC(=O)N1C(=O)NC(C)=CC1(C)C">(mol);
    validate_contains<"CC(C)CC(C(=O)NCC#N)C1=CC(=CC=C1)C2=CC=CC=C2">(mol);
    //validate_contains<"CC(C)CC(C)(C)[C@H](CCC#C)C(C)(C)C">(mol); // FIXME: stereo
    validate_contains<"CC(C)CC(CC(=O)NO)C(=O)NC(C(=O)NC(C)C(=O)NCCN)C(C)(C)C">(mol);
    validate_contains<"CC(C)CC(N)C#N">(mol);
    validate_contains<"CC(C)CC(N1C(=O)N2CCC3=C(NC4=C3C=CC=C4)C2(C)C1=O)C(=O)NC(CC(O)=O)C(O)=O">(mol);
    validate_contains<"CC(C)CC(NC(=O)C1=CC2=CC=CC=C2O1)C(=O)NC3CCOCC3=O">(mol);
    validate_contains<"CC(C)CC1=CC=C(CC2=CSCCS2)C=C1">(mol);
    validate_contains<"CC(C)CC1=CC=CC=C1">(mol);
    validate_contains<"CC(C)CCC1=CNC2=C1CCCC2">(mol);
    validate_contains<"CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C">(mol);
    //validate_contains<"CC(C)CN(C[C@@H](O)[C@H](CC1=CC=CC=C1)NC(=O)[C@@H]2CN(C(=O)O2)C3=CC=CC(=C3)C(F)(F)F)S(=O)(=O)C4=CC5=C(OCO5)C=C4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)C[C@H](NC(=O)N1CC(=O)NC2=CC=CC=C12)C(=O)N[C@@H](CC(O)=O)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(C)C[C@H](NC(=O)[C@@H](N)C(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](CCC(O)=O)C(N)=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)Cl">(mol);
    validate_contains<"CC(C)N(C(C)C)C1=CC=C(NC2=CC=NC3=CC(Cl)=CC=C23)C=C1O">(mol);
    validate_contains<"CC(C)N(C(O)C(C)(C)Cl)C1=C(C=CC=C1)C(C)(C)C(C)(C)C">(mol);
    validate_contains<"CC(C)N(C)C1=CC=CC2=C1C=CC=C2S(O)(O)NC3ONC=C3">(mol);
    validate_contains<"CC(C)N(CCC1=CNC2=C1C=CC=C2)C(C)C">(mol);
    //validate_contains<"CC(C)N1C(C=C[C@@H](O)C[C@@H](O)CC(O)=O)C(C2=CC=C(F)C=C2)C3=C1C=CC=C3">(mol); // FIXME: stereo
    validate_contains<"CC(C)N1N=C(C2=C1N=CN=C2N)C3=CC=CC=C3">(mol);
    validate_contains<"CC(C)NC(=O)NC=C">(mol);
    validate_contains<"CC(C)NC1=C(N=CC=C1)N2CCN(CC2)C(=O)C3=CC4=C(N3)C=CC(NS(C)(=O)=O)=C4">(mol);
    validate_contains<"CC(C)NCC(O)COC1=CC=CC2=C1C3=C(N2)C=CC=C3">(mol);
    validate_contains<"CC(C)NCC(O)COC1=CC=CC2=C1C=CN2">(mol);
    //validate_contains<"CC(C)NS(=O)(=O)C1=CC=C(N)C=C1.CCC(O)C(CC2=CC=CC=C2O)NC(=O)OC3CCOC3">(mol); // FIXME: components
    validate_contains<"CC(C)Nc1cccnc1N2CCN(CC2)C(=O)c3cc4cc(NS(C)(=O)=O)ccc4[nH]3">(mol);
    validate_contains<"CC(C)OC(=O)NCC1(CC(O)=O)CCCCC1">(mol);
    validate_contains<"CC(C)ON(CC(=O)NO)S(=O)(=O)OC1=CC=C(C=C1)C2=CC=C(C)C=C2">(mol);
    validate_contains<"CC(C)ON(CC(=O)NO)S(=O)(=O)SC1=CC=C(C=C1)C2=CC=C(C)C=C2">(mol);
    validate_contains<"CC(C)OP(=O)(OC(C)C)F">(mol);
    validate_contains<"CC(C)S(=O)(=O)N1C(N)=NC2=C1C=C(C=C2)C(=NO)C3=CC=CC=C3">(mol);
    //validate_contains<"CC(C)[C@@H](CO)N1C=C(C(O)=O)C(=O)c2cc(Cc3cccc(Cl)c3F)ccc12">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](N)C(=O)N1CCC[C@@H]1C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](C(C)C)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](N1CCCNC1=O)C(=O)N[C@H](C[C@H](O)[C@H](CC2=CC=CC=C2)NC(=O)COC3=C(C)C=CC=C3C)CC4=CC=CC=C4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)CNC(=O)[C@@H]1CCCN1)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)CON=C1/CC[C@]2(C)[C@@H]3CC[C@@]4(C)[C@@H](CC[C@]4(O)C#C)[C@H]3CCC2=C1)C(=O)N[C@@H](CC(O)=O)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)[C@H](C)NC(=O)OCc1ccccc1)C(=O)N[C@H](CO)Cc2ccccc2">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(O)=O)NC(C)=O)C(=O)N[C@@H](CC(O)=O)C=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)[N+](C)(C)CCC1=CNC2=C1C=CC=C2">(mol);
    validate_contains<"CC(C1CCC2C3CCC4CC(N)CCC4(C)C3CCC12C)N(C)C">(mol);
    validate_contains<"CC(C1NC2=CC=CC=C2NC1=O)[N+]([O-])=O">(mol);
    validate_contains<"CC(CC(O)=C)C(O)=O">(mol);
    validate_contains<"CC(CCC1C(C)C2(C)CCC1(C)O2)C3CCC(C)=CC3O">(mol);
    validate_contains<"CC(CO)C1=CC(CO)=CC=C1">(mol);
    validate_contains<"CC(Cc1ccccc1)NC=O">(mol);
    validate_contains<"CC(F)(F)C1=CC=CC=C1">(mol);
    //validate_contains<"CC(N)=C1/C(=O)NC2=CC=C(C)C=C12">(mol); // FIXME: stereo
    validate_contains<"CC(N)=C1C(=O)NC2=CC=C(C=C12)C3=CN=CO3">(mol);
}
