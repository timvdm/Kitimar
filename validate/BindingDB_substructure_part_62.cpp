#include "Validate.hpp"

void BindingDB_substructure_part_62(OpenBabel::OBMol &mol)
{
    // SMARTS 3051 - 3100
    validate_contains<"NCCC(=O)Nc1ccccc1">(mol);
    validate_contains<"NCCC1=CC(O)=C(O)C=C1">(mol);
    validate_contains<"NCCC1=CC=C(O)C=C1">(mol);
    validate_contains<"NCCCCC">(mol);
    validate_contains<"NCCCOC1=CC=CC=C1">(mol);
    validate_contains<"NCCNCCNCCN">(mol);
    validate_contains<"NCl">(mol);
    validate_contains<"NN1C=CC=C1C#N">(mol);
    validate_contains<"NNC1=CC=C(C=C1)C(O)=O">(mol);
    validate_contains<"NS(=O)(=O)c1ccc(cc1)C(=O)NCc2ccccc2">(mol);
    //validate_contains<"N[C@@H](CC1=CC2=C(N1)C=CC=C2)C(N)=O">(mol); // FIXME: stereo
    //validate_contains<"N[C@@H](CO)CC1=CC=CC=C1">(mol); // FIXME: stereo
    //validate_contains<"N[C@@H]([C@H]1CC[C@@H](CC1)NS(=O)(=O)c2ccc(F)cc2F)C(=O)N3CC[C@H](F)C3">(mol); // FIXME: stereo
    validate_contains<"Nc1c2CCCCc2nc3ccccc13">(mol);
    validate_contains<"Nc1nc(N)c(-c2cccc(Cl)c2)c(CCCc3ccccc3)n1">(mol);
    validate_contains<"Nc1nc(N)c2c(ccc3c4ccccc4ccc23)n1">(mol);
    validate_contains<"Nc1nc(N)c2nc(CN3c4ccccc4C=Cc5ccccc35)cnc2n1">(mol);
    validate_contains<"Nc1nc(N)c2nc(CNc3ccc4cc(Cl)c(Cl)cc4c3)cnc2n1">(mol);
    validate_contains<"Nc1ncccc1NCc2cc(ccc2OCc3ccccc3)-c4ccc5cc[nH]c5c4">(mol);
    validate_contains<"Nc1ncccn1">(mol);
    //validate_contains<"Nc1ncnc2n(cnc12)[C@@H]3O[C@H](COP(O)(=O)OP(O)(=O)NP(O)(O)=O)[C@@H](O)[C@H]3O">(mol); // FIXME: stereo
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CC3=CC4=CC=CC=C4C=C3OC2=O">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CC3=CC=CC=C3OC2=O">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CC=CC=C2NCC3=CC=NC=C3">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2CN(C3CCCCC3)C(=O)C2">(mol);
    validate_contains<"O=C(NN=CC1=CC=CC=C1)C2=NOC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"O=C(NN=CC1=CC=CC=C1)C2=NOC=C2">(mol);
    validate_contains<"O=C(NS(=O)(=O)C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CC2=C(N1)C=CC=C2">(mol);
    validate_contains<"O=C1CCCC(=O)N1">(mol);
    validate_contains<"O=C1CCCC2=CC=CC=C12">(mol);
    validate_contains<"O=C1CCCCN1">(mol);
    validate_contains<"O=C1CNC=C1">(mol);
    validate_contains<"O=C1CNCC2=CC=CC=C2N1">(mol);
    validate_contains<"O=C1CSC=N1">(mol);
    validate_contains<"O=C1Cc2c([nH]c3ccc(cc23)C#N)-c4ccccc4N1">(mol);
    validate_contains<"O=C1N=CN2N=CC=CC2=C1C3=CC=CC=C3">(mol);
    validate_contains<"O=C1NC(=O)C(=C1NC2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"O=C1NC2=C(C=CS2)C(=O)N1">(mol);
    validate_contains<"O=C1NC2=CC=CC=C2C1=O">(mol);
    validate_contains<"O=C1NC2=CC=CC=C2C=C1">(mol);
    validate_contains<"O=C1NCNC=C1">(mol);
    //validate_contains<"O=C1Nc2ccc(cc2C1=O)S(=O)(=O)N3CCC[C@H]3COc4cccnc4">(mol); // FIXME: stereo
    validate_contains<"O=C1OCCCO1">(mol);
    validate_contains<"O=C1OCCN1C2=CC=C(C=C2)N3CCOCC3=O">(mol);
    validate_contains<"O=C1c2ccccc2-c3n[nH]c4cccc1c34">(mol);
    validate_contains<"O=CC1=CNC2=CC=CC=C2C1=O">(mol);
    validate_contains<"O=CC1CNCCN1">(mol);
    validate_contains<"O=CN1CCC1">(mol);
    //validate_contains<"OC(=O)C1=C(C(O)=CC=C1)C(=O)C2=C(O)C=C(C=C2O)C(=O)O[C@H]3C4CC(C5CCCC45)[C@H]3NC(=O)C6=CC=C(O)C=C6">(mol); // FIXME: stereo
}
