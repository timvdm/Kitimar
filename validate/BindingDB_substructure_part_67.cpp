#include "Validate.hpp"

void BindingDB_substructure_part_67(OpenBabel::OBMol &mol)
{
    // SMARTS 3301 - 3350
    validate_contains<"O=C1CCNO1">(mol);
    validate_contains<"O=C1CCOC2=CC=CC=C12">(mol);
    validate_contains<"O=C1CNCCN1">(mol);
    validate_contains<"O=C1NC2=C(N1)C=CC=C2">(mol);
    validate_contains<"O=CNCCC1=CC=CC=C1">(mol);
    validate_contains<"O=S(=O)(NC1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"OC(=O)C1=CON=C1O">(mol);
    validate_contains<"OC(=O)CC(O)(CC(O)=O)C(O)=O">(mol);
    validate_contains<"OC(=O)CCC(=O)C(O)=O">(mol);
    validate_contains<"OC1=CC=CC=C1">(mol);
    validate_contains<"OCC1=CC2=C(C=CC=C2)N=C1Cl">(mol);
    //validate_contains<"O[C@@H]1[C@@H](O)[C@@H](CC2=CC=CC=C2)N(CC3=CC=CC=C3)C(=O)N(CC4=CC=CC=C4)[C@@H]1CC5=CC=CC=C5">(mol); // FIXME: stereo
    validate_contains<"S1C=CN=C1">(mol);
    validate_contains<"SC1=NC2=CC=CC=C2S1">(mol);
    validate_contains<"[H][N+]([H])([H])C">(mol);
    //validate_contains<"[NH3+][C@H](CC1=CC=CC=C1)C(=O)N2CCC[C@H]2C(=O)NCC3=CC=CC=C3">(mol); // FIXME: stereo
    validate_contains<"c2ccc1[nH]ccc1c2">(mol);
    validate_contains<"C1CCC2=C(C1)NC=C2">(mol);
    validate_contains<"CNC(N)=N">(mol);
    validate_contains<"COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC">(mol);
    validate_contains<"NC1=NC=CS1">(mol);
    validate_contains<"C(OC1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"C1(=NN=CS1)N">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=CC=NC=C2">(mol);
    validate_contains<"C1=CN=CC=N1">(mol);
    validate_contains<"C=CC">(mol);
    //validate_contains<"CC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc2ccc(cc2)[C@@H]3CC(=O)NS3(=O)=O)C(N)=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)(C)C1=CC=CC=C1">(mol);
    validate_contains<"CC(C)N1C=CC=C1">(mol);
    //validate_contains<"CC(C)[C@H](N1CCCNC1=O)C(=O)N[C@H](C[C@H](O)[C@H](Cc2ccccc2)NC(=O)COc3c(C)cccc3C)Cc4ccccc4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)N(C)CC1=CSC(=N1)C(C)C)C(=O)N[C@H](C[C@H](O)[C@H](CC2=CC=CC=C2)NC(=O)OCC3=CN=CS3)CC4=CC=CC=C4">(mol); // FIXME: stereo
    validate_contains<"CNC1=C2C=CC=CC2=NC=N1">(mol);
    validate_contains<"CNC1=CC=NC2=CC=CC=C12">(mol);
    validate_contains<"N#CC1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"NC(=N)C1=CC2=C(NC(=C2)C3=CC=CC=C3O)C=C1">(mol);
    validate_contains<"NC(N)=N">(mol);
    validate_contains<"NC1=CC=NC2=C1C=CC=C2">(mol);
    validate_contains<"NC1=NC(N)=NC(N)=N1">(mol);
    validate_contains<"NC1=NC(N)=NN1">(mol);
    validate_contains<"NCC1(CC(O)=O)CCCCC1">(mol);
    validate_contains<"O=C1CC2=CC=CC=C2N1">(mol);
    validate_contains<"O=C1CCCN1">(mol);
    validate_contains<"O=C1NC=NC2=CC=CC=C12">(mol);
    validate_contains<"OC(=O)CCCCC1CCSS1">(mol);
    //validate_contains<"[H]OC1=CC=C(C=C1)C([H])=C(/[H])C2=CC(O[H])=CC(O[H])=C2">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]12CC[C@]3(C)[C@@H](O)CC[C@@]3([H])[C@]1([H])CCC4=CC(O)=CC=C24">(mol); // FIXME: stereo
    validate_contains<"[O-][N+]">(mol);
    validate_contains<"C1=CC2=CC=CC=C2C=C1">(mol);
    validate_contains<"Br">(mol);
    validate_contains<"C1=CN2C=CC=CC2=N1">(mol);
}
