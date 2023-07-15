#include "Validate.hpp"

void BindingDB_substructure_part_37(OpenBabel::OBMol &mol)
{
    // SMARTS 1801 - 1850
    //validate_contains<"C[C@@H]([C@H](N)[C@@H](C=C)C1=CC=CC=C1)[C@@H]2CCCCC2">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H]1CN([C@@H](C)CN1C(=O)NC1=CC=C(F)C=C1F)C1=CC(=C(C=C1)C#N)C(F)(F)F">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H]1OC[C@H](N)[C@@H]1O">(mol); // FIXME: stereo
    //validate_contains<"C[C@@](COC1=CC=CC=C1Cl)(NC(=O)C2=CC=C(C=C2)C(F)(F)F)C#N">(mol); // FIXME: stereo
    //validate_contains<"C[C@H](C(=O)N1Cc2[nH]nc(NC(=O)c3ccc(cc3)N4CCN(C)CC4)c2C1)c5ccccc5">(mol); // FIXME: stereo
    //validate_contains<"C[C@H](CC(O)=O)NC(=O)[C@@H](N)CC1=CC=CC=C1">(mol); // FIXME: stereo
    //validate_contains<"C[C@H](N)C(=O)O[C@@H](C)C(=O)O[C@@H](C)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"C[C@H](OC1=CC=C2NC(=CC2=C1)C(=O)C3=CC4=C(N3)C=CC=C4)C5=CC=CC=C5">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1CN(C)C(=O)O1">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1CN(CCC(C(=O)NCC2=CC(=CC(=C2)C(F)(F)F)C(F)(F)F)C3=CSC(NC(=O)CC4=CC=CC=N4)=N3)CC[C@@]15C=CC6=CC=CC=C56">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1CNC(=O)O1">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1CNCC1C">(mol); // FIXME: stereo
    //validate_contains<"C[C@](O)(C1=CC=C(C=C1)C(=O)N(C2CC2)[C@@H]3CC[C@H](CCC#N)CC3)C(F)(F)F">(mol); // FIXME: stereo
    //validate_contains<"C[C@](O)(C1=CC=C(C=C1)C(=O)N(C2CC2)[C@H]3CC[C@@H](CC3)C4=CC=CN=C4)C(F)(F)F">(mol); // FIXME: stereo
    //validate_contains<"C[C@](O)(CSc1ccc(NC([CH2+])=O)cc1)C(N)=O">(mol); // FIXME: stereo
    validate_contains<"C[N+](C)(C)CC1COCO1">(mol);
    validate_contains<"C[N+](C)(CC=C)C1=CC=C(CCC(=C)CCC2=CC=C(C=C2)[N+](C)(C)CC=C)C=C1">(mol);
    validate_contains<"C[N+]1=CC=C(C=C1)C2=CC=[N+](C)C=C2">(mol);
    validate_contains<"C[N+]1=C[N+](C)=C[N+](C)=C1">(mol);
    validate_contains<"C[P+](C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"C[S](=O)=O">(mol);
    validate_contains<"C[Sb](C)C">(mol);
    validate_contains<"C[Se]C(N)N">(mol);
    //validate_contains<"C[Se]C(NC=O)=N/C=O">(mol); // FIXME: stereo
    //validate_contains<"C\C">(mol); // FIXME: stereo
    //validate_contains<"C\N=C\C1=CC=CC=C1O">(mol); // FIXME: stereo
    validate_contains<"Cc1cc2nc(N)[nH]c2cc1C">(mol);
    validate_contains<"Cc1ccc(cc1)S(=O)(=O)OCC(=O)Nc1ccc(C(O)=O)c(O)c1">(mol);
    validate_contains<"Cc1ccccc1-c2ccc(NC(=O)c3cc(nn3-c4ccc5onc(N)c5c4)C(F)(F)F)c(F)c2">(mol);
    validate_contains<"Cc1ccccc1N2C(=O)c3c(C)nc4ccc(Br)cc4c3C2=O">(mol);
    //validate_contains<"Cc1ccccc1NC(=O)Nc2ccc(CC(=O)N[C@@H](CCCCNC(=O)C=Cc3cccnc3)C(=O)N[C@@H](CCCC(O)=O)C(=O)NC4(CCCCC4)C(N)=O)cc2">(mol); // FIXME: stereo
    validate_contains<"Cc1csc(C(N)=O)c1Cl">(mol);
    //validate_contains<"Cc1n[nH]c2ccc(cc12)-c3cncc(OC[C@@H](N)Cc4cccc(OCC5CCCCC5)c4)c3">(mol); // FIXME: stereo
    validate_contains<"Cc1n[nH]cc1-c2ccccc2">(mol);
    validate_contains<"Cc1nc(O)c(O)c(n1)C(=O)NCc2ccc(F)cc2">(mol);
    validate_contains<"Cc2ccc1[nH]nc(N)c1c2">(mol);
    //validate_contains<"Cl.C1=CC=NC=C1">(mol); // FIXME: components
    //validate_contains<"Cl.CC1=C2C(N)C(=O)NC2=CC=C1">(mol); // FIXME: components
    //validate_contains<"Cl.Cl.C1=CC=CC=C1">(mol); // FIXME: components
    //validate_contains<"Cl.Cl.Cl.CN(CCCNC1=C2C=CC(Cl)=CC2=NC3=C1CCCC3)CCCNC4=C5C=CC(Cl)=CC5=NC6=C4C7CC(C6)C=C(C)C7">(mol); // FIXME: components
    validate_contains<"ClC1=CC=C(C=C1)C2(CCNCC2)C3=CC=C(C=C3)C4=CNN=C4">(mol);
    validate_contains<"ClC1=CC=C(C=C1)C2=CC=CO2">(mol);
    validate_contains<"ClC1=CC=C(C=C1)N2N=CC3=C2C(=O)NCC3">(mol);
    validate_contains<"ClC1=CC=C(C=C1)S(=O)(=O)NC2=CC(Cl)=CC(Cl)=C2">(mol);
    validate_contains<"ClC1=CC=C(NC(=O)C2=CC=CC=C2)C=C1">(mol);
    validate_contains<"ClC1=CC=CC(=C1)N2CCNCC2">(mol);
    validate_contains<"ClC1=CC=CC(Cl)=C1C(=O)NC2=CNN=C2C(=O)NC3CCNCC3">(mol);
    validate_contains<"ClC1=CC=NC(Cl)=N1">(mol);
    validate_contains<"ClC1=CN=CN=C1">(mol);
    validate_contains<"ClC1=NC2=C(C=CC=C2)C=C1C=O">(mol);
}
