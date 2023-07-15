#include "Validate.hpp"

void BindingDB_substructure_part_68(OpenBabel::OBMol &mol)
{
    // SMARTS 3351 - 3400
    validate_contains<"C1CC2CCOC2O1">(mol);
    validate_contains<"C1COC2=C(O1)C=CC=C2">(mol);
    validate_contains<"C1COC2=CC=CC=C2N1">(mol);
    //validate_contains<"CC(=O)N[C@@H]1[C@@H](O)C[C@@](O[C@H]1[C@H](O)[C@H](O)CO)(Oc2ccc3C(C)=CC(=O)Oc3c2)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"CC(O)=O">(mol);
    validate_contains<"CC1=CC=CC2=C1N=CC=C2">(mol);
    //validate_contains<"CCC1=CC2C[C@@H](C1)c3c(N)c4ccc(Cl)cc4nc3C2">(mol); // FIXME: stereo
    validate_contains<"CCCc1nc(C)c2C(=O)NC(=Nn12)c3cc(ccc3OCC)S(=O)(=O)N4CCN(CC)CC4">(mol);
    validate_contains<"COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4">(mol);
    validate_contains<"NC1=C2C=CC=CC2=CC=N1">(mol);
    //validate_contains<"NC1=CC=CC(CN2[C@H](CC3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H](CC4=CC=CC=C4)N(CC5=CC=CC(N)=C5)C2=O)=C1">(mol); // FIXME: stereo
    validate_contains<"NC1=CC=NN1">(mol);
    validate_contains<"NC1=NC2=CC=CC=C2N1">(mol);
    validate_contains<"NC1=NC=NC2=C1N=CN2">(mol);
    validate_contains<"NC1=NOC=C1">(mol);
    validate_contains<"O1C=NN=C1">(mol);
    //validate_contains<"OC1=CC=C(\C=C\C2=CC(O)=CC(O)=C2)C=C1">(mol); // FIXME: stereo
    validate_contains<"OC1=CC=NC2=C1C=CC=C2">(mol);
    validate_contains<"Oc1ccc(cc1)C=Cc2cc(O)cc(O)c2">(mol);
    validate_contains<"[O-][N+](=O)C1=CC=CC=C1">(mol);
    validate_contains<"c1ccccc1">(mol);
    validate_contains<"ClC1=CC=CC=C1">(mol);
    validate_contains<"C1=CC=C2C=C3C=CC=CC3=CC2=C1">(mol);
    //validate_contains<"C1=CC=CC=C1.C2=CC3=C(C=C2)C=CC=C3">(mol); // FIXME: components
    validate_contains<"C1CCNCC1">(mol);
    validate_contains<"C1CNC2=CC=CC=C2N1">(mol);
    validate_contains<"CC(=O)NO">(mol);
    //validate_contains<"CC1=C(NC=C1)C=C2/C(=O)NC3=C2C=CC=C3">(mol); // FIXME: stereo
    validate_contains<"CC1=CN=CS1">(mol);
    validate_contains<"CS(=O)c1ccc(cc1)-c2nc(-c3ccc(F)cc3)c([nH]2)-c4ccncc4">(mol);
    validate_contains<"N1C=CC2=CN=CN=C12">(mol);
    validate_contains<"NC1=CC=CC=C1">(mol);
    validate_contains<"NC1=CC=CC=N1">(mol);
    validate_contains<"NCC">(mol);
    validate_contains<"NS">(mol);
    validate_contains<"Nc1ccccc1">(mol);
    validate_contains<"O">(mol);
    validate_contains<"O1C=NC2=CC=CC=C12">(mol);
    validate_contains<"O=CNCCc1ccccc1">(mol);
    validate_contains<"C1CCCC1">(mol);
    //validate_contains<"CC(C)C[C@H](NC(=O)[C@@H]1CCCN1C(=O)[C@@H]([NH3+])C(C)C)C([O-])=O">(mol); // FIXME: stereo
    validate_contains<"O=C1OC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"OC(=O)C1=CC=CC=C1O">(mol);
    validate_contains<"C1CCC2CCCC2C1">(mol);
    validate_contains<"C1OC2=C(O1)C=CC=C2">(mol);
    validate_contains<"CCNC">(mol);
    validate_contains<"COc1cc2c(Nc3cc(CC(=O)Nc4cccc(F)c4F)[nH]n3)ncnc2cc1OCCCN(CCO)CC(C)C">(mol);
    validate_contains<"ClC1=CC=CC=C1Cl">(mol);
    validate_contains<"FC(F)(F)C1=CC=CC=C1">(mol);
    validate_contains<"N1C=CC2=C1N=CC=C2">(mol);
}
