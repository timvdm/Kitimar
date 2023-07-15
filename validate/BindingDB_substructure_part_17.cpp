#include "Validate.hpp"

void BindingDB_substructure_part_17(OpenBabel::OBMol &mol)
{
    // SMARTS 801 - 850
    validate_contains<"CC(=O)C(=O)CCC(O)=O">(mol);
    validate_contains<"CC(=O)C(C)=O">(mol);
    validate_contains<"CC(=O)C(CC(O)=O)NC(=O)C(CCCCNS(=O)(=O)C1=CC=C(O)C(=C1)C(O)=O)C2=CC=CC=C2">(mol);
    validate_contains<"CC(=O)C(CC1CCCC1)C2=CC=CC=C2">(mol);
    validate_contains<"CC(=O)C1(C)NCC(C)(C)O1">(mol);
    validate_contains<"CC(=O)C1=C(C=CC=C1)C(O)=O">(mol);
    validate_contains<"CC(=O)C1=CC(=CC(=C1)C(C)(C)C)C(C)(C)C">(mol);
    validate_contains<"CC(=O)C1=CC2=C(NC(C3CC=CC23)C4=CC5=C(OCO5)C=C4Br)C=C1">(mol);
    validate_contains<"CC(=O)C1=CC2=C(OCCO2)C=C1">(mol);
    validate_contains<"CC(=O)C1=CC=CC(=C1)C2=CN=C3C=CC(NCC4CC4)=NN23">(mol);
    validate_contains<"CC(=O)C1=CN=C(C)C=C1">(mol);
    validate_contains<"CC(=O)C1=CNC2=CC=CC=C2C1=O">(mol);
    validate_contains<"CC(=O)C=CC(=O)NO">(mol);
    validate_contains<"CC(=O)CC#N">(mol);
    validate_contains<"CC(=O)CC(=O)C(O)=O">(mol);
    validate_contains<"CC(=O)CC(C1=CC=CC=C1)C2=C(O)C3=CC=CC=C3OC2=O">(mol);
    validate_contains<"CC(=O)CC(NC(=O)C1=NC2=C(C=CC=C2)C=C1)C(=O)NC(CC3=CC=CC=C3)C(O)CN4CC5C(O)CCCC5CC4C(=O)NC(C)(C)C">(mol);
    validate_contains<"CC(=O)CC1=CC=CC=C1">(mol);
    validate_contains<"CC(=O)CC1CCOCC1">(mol);
    validate_contains<"CC(=O)CC=CC(=O)NO">(mol);
    validate_contains<"CC(=O)CCP(C)(=O)C1=CC=C(Br)C=C1">(mol);
    validate_contains<"CC(=O)N1CC(O)CC1C(O)=O">(mol);
    validate_contains<"CC(=O)N1CCN(C(C1)C(N)=O)C(C)=O">(mol);
    validate_contains<"CC(=O)N1CCN(CC1)CC2=CC(=CC=C2)C(=O)NC3=CC(NC4=NC(=CC=N4)C5=CC=CN=C5)=C(C)C=C3">(mol);
    validate_contains<"CC(=O)N1Cc2[nH]nc(NC(=O)c3ccc(C)cc3)c2C1">(mol);
    validate_contains<"CC(=O)NC(CC(O)=O)C(=O)NCC1=CC=CC=C1">(mol);
    validate_contains<"CC(=O)NC1(CC2=CC=CC=C2C1)NC(C)=O">(mol);
    validate_contains<"CC(=O)NC1=C(C=CC=C1)C(=O)NC2=CC=C(Cl)C=C2">(mol);
    validate_contains<"CC(=O)NC1=CC(C(C)=O)=C(N)C=C1">(mol);
    validate_contains<"CC(=O)NC1=CC(O)=C(C)C(C)=C1O">(mol);
    validate_contains<"CC(=O)NC1=CC=C(C)C=C1">(mol);
    validate_contains<"CC(=O)NC1=CC=C(C=C1)N2C=CN=C2">(mol);
    validate_contains<"CC(=O)NC1=CC=C(Cl)C=C1">(mol);
    validate_contains<"CC(=O)NC1=CC=C(NC2=CC(=O)C3=C(C=CC=C3)C2=O)C=C1">(mol);
    validate_contains<"CC(=O)NC1=CC=C(O)C=C1">(mol);
    validate_contains<"CC(=O)NC1=CC=CC=C1N">(mol);
    validate_contains<"CC(=O)NC1=NNC2=C1C=C(C=N2)C3=CC=CC=C3">(mol);
    validate_contains<"CC(=O)NC1=NNC2=C1C=NC=C2">(mol);
    //validate_contains<"CC(=O)NC1C(C=C(OC1C(O)C(O)CO)C(O)=O)N=C(/N)N">(mol); // FIXME: stereo
    validate_contains<"CC(=O)NC1COC2=C1SC3=CC=CC=C3N2">(mol);
    validate_contains<"CC(=O)NCC1CN(C(=O)O1)C2=CC=CC=C2">(mol);
    validate_contains<"CC(=O)NCCC1=CC=CC=C1">(mol);
    //validate_contains<"CC(=O)N[C@@H](CC(O)=O)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(=O)N[C@@H]1[C@@H](NC(N)=N)C=C(O[C@H]1[C@H](O)[C@H](O)CO)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CC(=O)Nc1ccc(SC[C@](C)(O)C(=O)NC2=CC=C(C(=C2)C(F)(F)F)N(=O)=O)cc1">(mol); // FIXME: stereo
    //validate_contains<"CC(=O)Nc1ccc(SC[C@](C)(O)C(=O)NC2CC=C(C(=C2)[C++](F)(F)(F)=C)N(=O)=O)cc1">(mol); // FIXME: stereo
    validate_contains<"CC(=O)OC1=CC=CC=C1">(mol);
    validate_contains<"CC(=O)OC1C2=C(C)C(CC(O)(C(OC(=O)C3=CC=CC=C3)C4C5(COC5CC(O)C4(C)C1=O)OC(C)=O)C2(C)C)OC(=O)C(O)C(NC(=O)C6=CC=CC=C6)C7=CC=CC=C7">(mol);
    validate_contains<"CC(=O)OC1CC2C3CCC1C23">(mol);
    validate_contains<"CC(=O)OC1CC2C3OC3C(C1)[N+]2(C)C">(mol);
}
