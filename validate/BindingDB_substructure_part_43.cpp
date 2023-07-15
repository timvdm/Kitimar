#include "Validate.hpp"

void BindingDB_substructure_part_43(OpenBabel::OBMol &mol)
{
    // SMARTS 2101 - 2150
    validate_contains<"NC1=NC=NC2=C1C3=CC=CC=C3N2">(mol);
    validate_contains<"NC1=NC=NC2=C1N=CN2C3OC(CO)C(OP(O)(O)=O)C3O">(mol);
    //validate_contains<"NC1=NC=NC2=C1N=CN2[C@@H]3O[C@@H](C[C@H]3O)OP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O">(mol); // FIXME: stereo
    validate_contains<"NC1=NC=NN1C=O">(mol);
    validate_contains<"NC1=NN(C(=O)C2=CC=CC=C2)C(N)=N1">(mol);
    validate_contains<"NC1=NN(C=O)C(N)=N1">(mol);
    validate_contains<"NC1=NNC2=C1C=CS2">(mol);
    validate_contains<"NC1=NSC=C1">(mol);
    validate_contains<"NC1=NSN=C1N">(mol);
    validate_contains<"NC1C2CC3CC(C2)CC1C3">(mol);
    validate_contains<"NC1CC1">(mol);
    validate_contains<"NC1CC2=CC=CC3=C2C4C(CC=CC14)O3">(mol);
    validate_contains<"NC1CCC=C(C1)C(O)=O">(mol);
    validate_contains<"NC1CCC=CC1">(mol);
    validate_contains<"NC1CCN(CC1)C(=O)C2=CC=C(NC3=NC=CC(=N3)C4=CC5=C(S4)C=CC=C5)C=C2">(mol);
    validate_contains<"NC1N=C2C=CC=CC2=N1">(mol);
    validate_contains<"NC1N=CC=S1">(mol);
    validate_contains<"NC1NC2=CN=CN=C2N1">(mol);
    validate_contains<"NC=CCC(P(O)(O)=O)P(O)(O)=O">(mol);
    validate_contains<"NC=N">(mol);
    validate_contains<"NCC(=O)NCC(=O)NCC(O)=O">(mol);
    validate_contains<"NCC(C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    //validate_contains<"NCC(O)=O.C[C@H](N)C(O)=O.OC(=O)[C@@H]1CCCN1">(mol); // FIXME: components
    validate_contains<"NCC1(CCCCC1)CC(O)=O">(mol);
    validate_contains<"NCC1=CC(=O)C2=CC=CC=C2O1">(mol);
    validate_contains<"NCC1=CC2=C(CCCC2)C=C1">(mol);
    validate_contains<"NCC1=CC2=C(CCCN2)C=C1">(mol);
    validate_contains<"NCC1=CC=C(CN2CCOCC2)C=C1">(mol);
    validate_contains<"NCC1=CC=CC(=C1)C(O)=O">(mol);
    validate_contains<"NCC1=CC=CC=C1C(O)=O">(mol);
    validate_contains<"NCC1CC1">(mol);
    validate_contains<"NCC1CN(C(=O)O1)C2=CC=CC=C2">(mol);
    validate_contains<"NCC1COCO1">(mol);
    validate_contains<"NCCC(O)=O">(mol);
    validate_contains<"NCCC1=CN(C2=C1C=CC=C2)C3=CC=CC=C3">(mol);
    //validate_contains<"NCCC1CCN(CC1)C(=O)[C@H](Cc2cccc(c2)C(N)=[NH2+])NS(=O)(=O)c3cccc4ccccc34">(mol); // FIXME: stereo
    validate_contains<"NCCCCC(N)C(=O)N1CC(CC1C#N)N=[N+]=[N-]">(mol);
    validate_contains<"NCCCCC(N)C(=O)N1CCCC1C#N">(mol);
    validate_contains<"NCCN1C(=O)C2=CC(=CC3=C2C(=CC(=C3N)S(O)(=O)=O)C1=O)S(O)(=O)=O">(mol);
    validate_contains<"NCCc1c(F)cccc1Cl">(mol);
    validate_contains<"NCCc1c[nH]c2ccc(O)cc12">(mol);
    validate_contains<"NCCc1ccccc1">(mol);
    validate_contains<"NCNC1CCCCN(CC(=O)N2CCCC2)C1=O">(mol);
    //validate_contains<"NC[C@@H](NC(=O)C1=CC=C(S1)C2=CC=NC3=C2C=CN3)C4=CC=CC=C4">(mol); // FIXME: stereo
    //validate_contains<"NC[C@H](NC(=O)C1=CC=C(S1)C2=CC=NC3=C2C=CN3)C4=CC=CC=C4">(mol); // FIXME: stereo
    validate_contains<"NN(CC(=O)N1CCSC1)C2CCN(CC(=O)NC3=NC=CC=C3)CC2">(mol);
    validate_contains<"NN1CC(=O)NC1=O">(mol);
    validate_contains<"NNC1=NN=CC2=C1C=CC=C2">(mol);
    validate_contains<"NOOCC1=CC(COON)=CC=C1">(mol);
    validate_contains<"NS(=O)(=O)C1=CC(=CC=C1O)C(O)=O">(mol);
}
