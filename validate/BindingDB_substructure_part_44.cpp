#include "Validate.hpp"

void BindingDB_substructure_part_44(OpenBabel::OBMol &mol)
{
    // SMARTS 2151 - 2200
    validate_contains<"NS(=O)(=O)C1=CC=C(C=C1)C2=CNC3=C2C=CC(NC(=O)C4CCNCC4)=N3">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=C(C=C1)C2=CSC=N2">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=C(Cl)C=C1">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=C(NC(=O)CSC2=NN=C(Br)N2C3=C4C=CC=CC4=C(C=C3)C5CC5)C(Cl)=C1">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=C(NC2=NC=CC=N2)C=C1">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=C(NS(=O)(=O)C2=CC=CC=C2)S1">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=CC2=C1C=CN=C2">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=CC=C1C2=CC=C(CC(=O)C3=CC=NN3C4=CC=C5C=CC=CC5=C4)C=C2">(mol);
    validate_contains<"NS(=O)(=O)C1=CC=CC=C1C2=CC=CC=C2">(mol);
    validate_contains<"NS(=O)(=O)C1=NC=CC=C1">(mol);
    validate_contains<"NS(=O)(=O)C1=NN=C(NS(=O)(=O)C2=CC=CC=C2)S1">(mol);
    validate_contains<"NS(=O)(=O)C1=NN=C(S1)C2=C(Cl)C=CC=C2">(mol);
    validate_contains<"NS(=O)(=O)Nc1ccc2N=C(NS(=O)(=O)c2c1)c1c(O)c2ccccc2n(NC2CCC2)c1=O">(mol);
    validate_contains<"NS(=O)(=O)c1ccc(cc1)C(=O)NNC(=O)Nc2ccc(Cl)c(Cl)c2">(mol);
    validate_contains<"NS(=O)(=O)c1cccc(Nc2nccc(Nc3c(Cl)ccc4OCOc34)n2)c1">(mol);
    validate_contains<"NS(=O)=O">(mol);
    //validate_contains<"N[C@@H](CC(=O)N1CCN2C(C1)=NN=C2C(F)(F)F)CC3=C(F)C=C(F)C(F)=C3">(mol); // FIXME: stereo
    //validate_contains<"N[C@@H](CC1=CC=CC=C1)[C@H](O)C(=O)N2CCSC2">(mol); // FIXME: stereo
    //validate_contains<"N[C@@H](CO)C=CSC1=CC=C2C=CC=CC2=C1">(mol); // FIXME: stereo
    //validate_contains<"N[C@@H](Cc1ccc(cc1)-c2ccc(cc2)C#N)C(=O)N3CCSC3">(mol); // FIXME: stereo
    //validate_contains<"N[C@H](CC1=CC2=C(N1)C=CC=C2)C(N)=O">(mol); // FIXME: stereo
    //validate_contains<"N[C@H](COC1=CN=CC(=C1)C2=CC3=C(NC(=O)C3=C/C4=CC=CO4)C=C2)CC5=CNC6=C5C=CC=C6">(mol); // FIXME: stereo
    //validate_contains<"N[C@H](COC1=CN=CC(=C1)C2=CC=C3C=NC=CC3=C2)CC4=CNC5=C4C=CC=C5">(mol); // FIXME: stereo
    //validate_contains<"N[C@H]1CC=CCC1C2=CC=CC=C2Cl">(mol); // FIXME: stereo
    //validate_contains<"N[C@H]1CC=CC[C@H]1C2=CC=CC=C2Cl">(mol); // FIXME: stereo
    //validate_contains<"N[C@H]1COC(=C)[C@H]1O">(mol); // FIXME: stereo
    validate_contains<"Nc1c2CCCCc2nc2ccccc12">(mol);
    validate_contains<"Nc1ccc(P)cc1">(mol);
    validate_contains<"Nc1ccc(cc1)S(=O)(=O)N2CCCC2">(mol);
    validate_contains<"Nc1ccc2c(Nc3cccc(Br)c3)ncnc2c1">(mol);
    validate_contains<"Nc1nc(N)c2nccnc2n1">(mol);
    validate_contains<"Nc1nc(Nc2ccc(cc2)S(N)(=O)=O)nn1C(=S)Nc3c(F)cccc3F">(mol);
    validate_contains<"Nc1nccc(n1)-c2cc(Cl)sc2Cl">(mol);
    //validate_contains<"Nc1ncnc2n(cnc12)[C@@H]3O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]3O">(mol); // FIXME: stereo
    validate_contains<"Nc1ncnc2n(nc(-c3ccccc3)c12)-c4ccccc4">(mol);
    //validate_contains<"O.CC(=O)NC1=CC=C(NC(C)=O)C=C1">(mol); // FIXME: components
    //validate_contains<"O.CC12CCC3C(CCC4=CC(O)=CC=C34)C1CCC2O">(mol); // FIXME: components
    validate_contains<"O1C(=CC=C1C2=CC=CC=C2)N3C=CC=C3">(mol);
    validate_contains<"O1C2=C(C=CC=C2)N=C1C1=CC=CC=C1">(mol);
    validate_contains<"O1C2=C(OC3=C1C=CC=C3)C=CC=C2">(mol);
    validate_contains<"O1C2=CC=CC=C2C3=C1C=CC=C3">(mol);
    validate_contains<"O1C2=CC=CC=C2N=C1C3=CC=CC=C3">(mol);
    validate_contains<"O1C2=CC=CC=C2N=CC3=C1C=CC=C3">(mol);
    validate_contains<"O1C=CC2=C1C=CC=C2">(mol);
    validate_contains<"O1C=CC2=C1OC3=C2C=CO3">(mol);
    validate_contains<"O1C=CC2=CC3=C(C=CO3)C=C12">(mol);
    validate_contains<"O1C=NC2=CN=CC=C12">(mol);
    validate_contains<"O1N=CC2=CC=CN=C12">(mol);
    validate_contains<"O1c2ccc(cc2(OC1))CC(COC)C(COC)Cc3ccc4OCOc4(c3)">(mol);
    validate_contains<"O=C(C(=O)C1=CC=CC=C1)C2=CC=CC=C2">(mol);
}
