#include "Validate.hpp"

void BindingDB_substructure_part_51(OpenBabel::OBMol &mol)
{
    // SMARTS 2501 - 2550
    //validate_contains<"OCC1CC(O)C(O)CO1.OCC2CC(O)C(O)CO2.OCC3CC(O)C(O)CO3.OCC4CC(O)C(O)CO4.OCC5CC(O)C(O)CO5.OCC6CC(O)C(O)CO6.OCC7CC(O)C(O)CO7">(mol); // FIXME: components
    validate_contains<"OCC1CCC(O)C1O">(mol);
    validate_contains<"OCCC1=CC=CC=C1">(mol);
    validate_contains<"OCCC1=CNC2=C1C=CC=C2">(mol);
    validate_contains<"OCCCNC1=NC=CC(=C1)C1=NC(NC2=CC=CC(Cl)=C2)=NC=C1">(mol);
    validate_contains<"OCCN(CCO)c1-nc(C2CCCCC2)-c3-nc(-nc(C4CCCCC4)-c3-n1)N(CCO)CCO">(mol);
    validate_contains<"OCCN1C=C(C(O)=O)C(=O)c2cc(Cc3cccc(Cl)c3F)ccc12">(mol);
    validate_contains<"OCCNC(=O)C1=CC2=NN3C(OCC4=CC=CC=C34)=C2C=C1">(mol);
    //validate_contains<"OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](O)[C@H](O[C@@H]4CO)O[C@H]5[C@H](O)[C@@H](O)[C@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](O)[C@H](O[C@@H]6CO)O[C@H]7[C@H](O)[C@@H](O)[C@H](O[C@@H]7CO)O[C@H]1[C@H](O)[C@H]2O">(mol); // FIXME: stereo
    validate_contains<"ON(O)C1=CC2=C(C=NN2)C=C1">(mol);
    validate_contains<"ON(O)C1=CC=C2CN=NC2=C1">(mol);
    validate_contains<"ON=C1C(=C)NC2=C1C=CC=C2">(mol);
    //validate_contains<"ON=C1C(NC2=C1C=CC=C2)=C3/C=CNC3=O">(mol); // FIXME: stereo
    validate_contains<"ONC(=O)CCCCCCC(=O)Nc1cccc(c1)-c1ccccc1">(mol);
    validate_contains<"ONC(=O)CCCCCN1C(=O)C2=CC=CC3=C2C(=CC=C3)C1=O">(mol);
    validate_contains<"ONC(=O)CCCP(O)(O)=O">(mol);
    validate_contains<"OP(O)(=O)C(CC=C[N+]1=CC=CC=C1)P(O)(O)=O">(mol);
    validate_contains<"OP(O)(=O)C(F)(F)C1=CC=CC=C1">(mol);
    validate_contains<"OP(O)(=O)OC1=CC=CC=C1">(mol);
    validate_contains<"OS">(mol);
    validate_contains<"OS(=O)(=O)C1=CC=CC=C1">(mol);
    validate_contains<"OS(O)(=O)=O">(mol);
    validate_contains<"OS(O)(O)O">(mol);
    //validate_contains<"O[C@@H]1[C@@H](O)[C@@H](CC2=CC=CC=C2)N(CC3=CC(=CC=C3)C(=O)NC4=CC=CC=N4)C(=O)N(CC5=CC=CC(=C5)C(=O)NC6=NC=CC=C6)[C@@H]1CC7=CC=CC=C7">(mol); // FIXME: stereo
    //validate_contains<"O[C@@H]1[C@@H](O)[C@@H](CC2=CC=CC=C2)N(CC3=CC=C4C=CC=CC4=C3)C(=O)N(CC5=CC6=CC=CC=C6C=C5)[C@@H]1CC7=CC=CC=C7">(mol); // FIXME: stereo
    //validate_contains<"O[C@H]1CCCC[C@H]1O">(mol); // FIXME: stereo
    //validate_contains<"O[C@H]1CN2CCC1CC2">(mol); // FIXME: stereo
    //validate_contains<"O[C@H]1[C@@H](O)[C@H](O[C@@H]1CO[P@](O)(=O)OP(O)(O)=O)[N+]2=CCC(=O)NC2=O">(mol); // FIXME: stereo
    validate_contains<"O[N+](=O)C1=CC2=C(CN=C2)C=C1">(mol);
    validate_contains<"Oc1c(O)c(O)c2c(=O)cc(oc2c1)c1ccccc1">(mol);
    validate_contains<"Oc1c(nc(N2CCCCS2(=O)=O)c3cccnc13)-c4nnc(Cc5ccc(F)cc5)[nH]4">(mol);
    validate_contains<"Oc1cc(O)c2C(=O)C(O)=C(Oc2c1)c3ccc(O)c(O)c3">(mol);
    validate_contains<"Oc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O">(mol);
    validate_contains<"Oc1nc(Cc2ccccc2)nc(C(=O)NCc3ccc(F)cc3)c1O">(mol);
    validate_contains<"Oc1nc(nc(C(=O)NCc2ccc(F)cc2)c1O)-c3cccs3">(mol);
    validate_contains<"Oc1ncnc(C(=O)NCc2ccc(F)cc2)c1O">(mol);
    validate_contains<"PC1=CC=CC=C1">(mol);
    validate_contains<"PCCNC1=C(CC=C)NC=C1">(mol);
    validate_contains<"S">(mol);
    //validate_contains<"S.N1C=CC=C1.C2=CC3=C(C=C2)C=CC=C3">(mol); // FIXME: components
    validate_contains<"S1C(=CN=C1C2=CC=NC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"S1C2=C(C=CC=C2)N=C1C3=CC=CC=C3">(mol);
    validate_contains<"S1C=CC2=C1C=CC=C2">(mol);
    validate_contains<"S1C=CN2C=C(N=C12)C1=CC=CC=C1">(mol);
    validate_contains<"S1C=NC(=C1)N1C=CC=C1">(mol);
    validate_contains<"S1C=NC(=C1)N1C=CN=C1">(mol);
    validate_contains<"S1C=NC2=C1C=NC=N2">(mol);
    validate_contains<"S1C=NN=C1">(mol);
    validate_contains<"S1N=C2C=CC=CC2=N1">(mol);
    validate_contains<"S=C(NCCC1=CC=CC=C1)NC2=NC=CC=C2">(mol);
}
