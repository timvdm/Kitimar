#include "Validate.hpp"

void BindingDB_substructure_part_65(OpenBabel::OBMol &mol)
{
    // SMARTS 3201 - 3250
    //validate_contains<"C[C@@]12CCCC1C3CCC4CCCCC4C3CC2">(mol); // FIXME: stereo
    validate_contains<"Cc1nc(Nc2ncc(s2)C(=O)Nc3c(C)cccc3Cl)cc(n1)N4CCN(CCO)CC4">(mol);
    validate_contains<"Cl">(mol);
    validate_contains<"ClC1=CC=C(NC=O)C=C1">(mol);
    validate_contains<"ClC1=CC=CC(OCC2=CC=CC=C2)=C1Cl">(mol);
    validate_contains<"FC(F)(F)c1ccccc1S(=O)(=O)NNC(=O)c2cc3cc(Cl)ccc3o2">(mol);
    validate_contains<"FC1=CC=CC=C1">(mol);
    validate_contains<"N1C=CC2=C1C=C3C=CC=CC3=C2">(mol);
    validate_contains<"N1C=CC2=CC=CN=C12">(mol);
    //validate_contains<"N1C=CC=C1.C2CCCCC2">(mol); // FIXME: components
    validate_contains<"N1C=NC2=C1C=NC=N2">(mol);
    validate_contains<"N1C=NC2=CN=CN=C12">(mol);
    validate_contains<"N1C=NN=C1">(mol);
    validate_contains<"N1N=CC2=C1N=CN=C2">(mol);
    validate_contains<"NC(=N)C1=CC2=C(NC(=N2)C3=CC=CC=C3O)C=C1">(mol);
    validate_contains<"NC(=N)C1=CC=CC=C1">(mol);
    //validate_contains<"NC(=[NH2+])C1=CC=C(CNC(=O)[C@@H]2CCCN2C(=O)[C@H]([NH3+])C(C3=CC=CC=C3)C4=CC=CC=C4)C=C1">(mol); // FIXME: stereo
    validate_contains<"NC1=NC2=C(C=CC=C2)C(N)=N1">(mol);
    validate_contains<"NC1=NC=NC2=C1C(C3=CC=CC=C3)=NN2C(C)(C)C">(mol);
    validate_contains<"NCC1(CC1)C(O)=O">(mol);
    validate_contains<"NCCC1=CC=CC=C1">(mol);
    validate_contains<"NCCC1=CN=CN1">(mol);
    validate_contains<"NNC1=CC=CC=C1">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=CC=C1">(mol);
    //validate_contains<"N[C@@H](CC(=O)N1CCn2c(C1)nnc2C(F)(F)F)Cc3cc(F)c(F)cc3F">(mol); // FIXME: stereo
    validate_contains<"O1C=CN=C1">(mol);
    validate_contains<"O=C(O)CCc1cc(O)cc(O)c1">(mol);
    validate_contains<"O=C1CN=CC2=CC=CC=C2N1">(mol);
    validate_contains<"O=C1COC(=O)N1CC2=CC=CC=C2">(mol);
    validate_contains<"O=C1NC=CC=C1">(mol);
    validate_contains<"O=C1NCCC2=C1C=CN2">(mol);
    validate_contains<"O=N(=O)C1=CC=CC=C1">(mol);
    validate_contains<"OC(=O)C=O">(mol);
    validate_contains<"OC1=CC(=O)OCC1">(mol);
    validate_contains<"ON=O">(mol);
    validate_contains<"OS(=O)(=O)Nc1cccc(CCC(=O)N2CCc3ccc(NS(O)(=O)=O)cc3C2)c1">(mol);
    //validate_contains<"S.N1C=CC=C1.C2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"S1C=CC2=CC=CC=C12">(mol);
    validate_contains<"S1C=CC=C1">(mol);
    validate_contains<"S1C=CC=C1N1C=NC2=CC=CC=C12">(mol);
    validate_contains<"SC1=NC=CC=N1">(mol);
    validate_contains<"[H]C1=CN=CN=C1">(mol);
    validate_contains<"[H]C1=NC2=CC=CC=C2C(NC3=CC=CC=C3)=C1[H]">(mol);
    validate_contains<"[H]N([H])S(C)(=O)=O">(mol);
    validate_contains<"[H]OC(=O)CCNC(C)=O">(mol);
    //validate_contains<"[H]O[C@@]1([H])[C@@]([H])(O[C@]([H])(C([H])([H])O[P@@]([O-])(=O)O[P@@]([O-])(=O)OC([H])([H])C([H])([H])C2=C(N=C(S2)C([O-])=O)C([H])([H])[H])[C@@]1([H])O[H])N1C([H])=NC2=C1N=C([H])N=C2N([H])[H]">(mol); // FIXME: stereo
    //validate_contains<"[H][C@](CC1=CC=C(C=C1)C(N)[NH3+])(NC(=O)CNS(=O)(=O)C2=CC=C3C=CC=CC3=C2)C=O">(mol); // FIXME: stereo
    //validate_contains<"[NH3+][C@H](C(C1=CC=CC=C1)C2=CC=CC=C2)C(=O)N3CCC[C@H]3C(=O)NCC4=CC=CC(Cl)=C4">(mol); // FIXME: stereo
    validate_contains<"[O-]C(=O)CCCCC1CCSS1">(mol);
    validate_contains<"[O-][NH+]=O">(mol);
}
