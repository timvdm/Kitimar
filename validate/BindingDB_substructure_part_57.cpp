#include "Validate.hpp"

void BindingDB_substructure_part_57(OpenBabel::OBMol &mol)
{
    // SMARTS 2801 - 2850
    validate_contains<"C1CC2=CC=CC=C2SC2=CC=CC=C12">(mol);
    validate_contains<"C1CC2CCC3C(CCC4CCCCC34)C2C1">(mol);
    validate_contains<"C1CC2CCCC3=C2C(C1)=CC=C3">(mol);
    validate_contains<"C1CCC(CC1)COC2=CC=CC=C2">(mol);
    validate_contains<"C1CCC2(C1)CCC3(C2)CCCCC3">(mol);
    validate_contains<"C1CCC2(C1)CCCC2">(mol);
    validate_contains<"C1CCC2(CC1)CCCCC2">(mol);
    validate_contains<"C1CCC2=C(C1)C=CC3=C2C=CC=C3">(mol);
    validate_contains<"C1CCC2=C(C1)C=CC=C2">(mol);
    validate_contains<"C1CCC2C(C1)CCC3C4CCCC4CCC23">(mol);
    validate_contains<"C1CCC2C(C1)CCC3CCCCC23">(mol);
    validate_contains<"C1CCC2CC3CCCCC3CC2C1">(mol);
    validate_contains<"C1CCC2CCCCC2C1">(mol);
    //validate_contains<"C1CCCC1.C2=CC3=C(C=C2)C=CC=C3">(mol); // FIXME: components
    //validate_contains<"C1CCCC1.N2C=CC=C2">(mol); // FIXME: components
    //validate_contains<"C1CCCCC1.C2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"C1CCOCC1">(mol);
    validate_contains<"C1CN=C(N1)C2=CC=CC=C2">(mol);
    validate_contains<"C1CNC2=CC=CC=C2C1">(mol);
    validate_contains<"C1CNCNC1">(mol);
    validate_contains<"C1C[N:1]CNC1">(mol);
    validate_contains<"C1NCC(=C1)C2CNC3=CC=CC=C23">(mol);
    validate_contains<"C1NCC2=C1C3=C(NC4=CC=CC=C34)C5=C2C6=C(N5)C=CC=C6">(mol);
    validate_contains<"C1NCC2=C1C=CC=C2">(mol);
    validate_contains<"C1NCC2=CC=CC=C12">(mol);
    validate_contains<"C1NN=C2C1CC3=CC=CC=C23">(mol);
    validate_contains<"C1NNNN1">(mol);
    validate_contains<"C1NOC=C1">(mol);
    validate_contains<"C=C1C(=O)NC2=CC=CC=C12">(mol);
    validate_contains<"C=CC(=O)NC1=CC2=C(C=C1)N=CN=C2NC3=CC=CC=C3">(mol);
    validate_contains<"C=CN1C=CCC1=O">(mol);
    validate_contains<"CC#N">(mol);
    //validate_contains<"CC(=NN)C(C)=N/O">(mol); // FIXME: stereo
    validate_contains<"CC(=O)C1=C(O)C=CC=C1">(mol);
    validate_contains<"CC(=O)N1Cc2[nH]nc(NC(=O)c3ccc(F)cc3)c2C1">(mol);
    validate_contains<"CC(=O)N1Cc2[nH]nc(NC(=O)c3ccc(cc3)C#N)c2C1">(mol);
    validate_contains<"CC(=O)Nc1ccc(cc1)C(=O)Nc2ccccc2N">(mol);
    //validate_contains<"CC(C)(C)NC(=O)[C@@H]1CC2CCCCC2CN1C[C@@](O)(CC3=CC=CC=C3)NC(=O)[C@H](CC(N)=O)NC(=O)C4=CC=C5C=CC=CC5=N4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)(C)NC(=O)[C@@H]1CCCCN1C[C@@H](O)[C@@H]2Cc3ccc(OCCCC(=O)N[C@@H](CC(N)=O)C(=O)N2)cc3">(mol); // FIXME: stereo
    //validate_contains<"CC(C)(C)NC(=O)[C@@H]1C[C@@H]2CCCC[C@@H]2CN1C[C@@H](O)[C@H](Cc3ccccc3)NC(=O)[C@H](CC(N)=O)NC(=O)c4ccc5ccccc5n4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)(C)OC(=O)N[C@@H](CCCCCS)C(=O)NC1CCCCCC1">(mol); // FIXME: stereo
    validate_contains<"CC(C)(CO)CC1=CC(=C(O)C(=C1)C(C)(C)C)C(C)(C)C">(mol);
    //validate_contains<"CC(C)(Oc1ccccc1)C(=O)N[C@@H]2[C@@H]3C[C@H]4CC2C[C@@](C4)(C3)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)C">(mol);
    //validate_contains<"CC(C)C1(C)SC(N[C@@H](C)c2ccccc2Br)=NC1=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)C1=CC=CC=C1">(mol);
    validate_contains<"CC(C)CC1=CC=C(C=C1)C(C)C(O)=O">(mol);
    //validate_contains<"CC(C)CCCCCN1[C@H](Cc2ccccc2)[C@H](O)[C@@H](O)[C@@H](Cc3ccccc3)N(CCCCCC(C)C)C1=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)CCN1C(C)C(=O)N(C)c2cnc(Nc3cc(F)c(O)c(F)c3)nc12">(mol);
    //validate_contains<"CC(C)CN(C[C@@H](O)[C@H](Cc1ccccc1)NC(=O)O[C@H]2CCOC2)S(=O)(=O)c3ccc(N)cc3">(mol); // FIXME: stereo
}
