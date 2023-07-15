#include "Validate.hpp"

void BindingDB_substructure_part_21(OpenBabel::OBMol &mol)
{
    // SMARTS 1001 - 1050
    validate_contains<"CC1(C)CCCC2=C1C=CC(C=O)=C2">(mol);
    validate_contains<"CC1(C)OC2COC3(CN)OC(C)(C)OC3C2O1">(mol);
    validate_contains<"CC1(C)SCNC1C(N)=O">(mol);
    //validate_contains<"CC1(CCN(CC1)C2=NC=CC(=C2)C(O)=O)NCC(=O)N3[C@H](CC[C@H]3C#N)C#C">(mol); // FIXME: stereo
    validate_contains<"CC1(COC(OC1)C2=NC(C3=CC=NC(=N3)C4CC4)=C(N2)C5=CC=C(F)C=C5)C(=O)N6CCOCC6">(mol);
    validate_contains<"CC1(COC(OC1)C2=NC(C3=CC=NC(N)=N3)=C(N2)C4=CC=C(F)C=C4)C(=O)N5CCOCC5">(mol);
    validate_contains<"CC1(COC(OC1)C2=NC(C3=CC=NC(NC4CC4)=N3)=C(N2)C5=CC=C(F)C=C5)C(=O)N6CCOCC6">(mol);
    validate_contains<"CC1(COC(OC1)C2=NC(C3=CC=NC=C3)=C(N2)C4=CC=C(F)C=C4)C(=O)N5CCOCC5">(mol);
    validate_contains<"CC1(COC(OC1)C2=NC(C3=CC=NC=N3)=C(N2)C4=CC=C(F)C=C4)C(=O)N5CCOCC5">(mol);
    validate_contains<"CC12CC(=O)C3C(CCC4=CC(=O)CCC34C)C1CCC2(O)C(=O)CO">(mol);
    validate_contains<"CC12CC(O)CC1C3CCC4=CC(O)=CC=C4C3CC2">(mol);
    validate_contains<"CC12CC1(C)CN(C2)C3=CC=CC=C3">(mol);
    validate_contains<"CC12CC3CCC4=CC(O)=CC=C4C3CC1CCC2O">(mol);
    validate_contains<"CC12CCC(=O)C=C1CCC3C4CCC(C(=O)CO)C4(CO)CC(O)C23">(mol);
    validate_contains<"CC12CCC3C(CCC4=CC(O)=CC=C34)C1CC(O)C2O">(mol);
    validate_contains<"CC12CCC3C(CCC4=CC(O)=CC=C34)C1CCC2=O">(mol);
    validate_contains<"CC12CCC3C(CCC4=CC=CC=C34)C1CCC2O">(mol);
    validate_contains<"CC12CCC3C(CCC4=CC=CC=C34)C1CCCC2O">(mol);
    validate_contains<"CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O">(mol);
    validate_contains<"CC12CCCC1C3CCC4=CC(O)=CC=C4C3CC2O">(mol);
    validate_contains<"CC12CCCC1C3CCC4=CC=CC=C4C3CC2">(mol);
    validate_contains<"CC12OC(CC1(O)C(O)=O)N3C4=C(C=CC=C4)C5=C3C6=C(C7=C5C(=O)NC7)C8=C(C=CC=C8)N26">(mol);
    validate_contains<"CC1=C(C#N)C2=C(C=CC=C2)C=N1">(mol);
    //validate_contains<"CC1=C(C(=O)NCC2=CC=NC=C2)C(C)=C(N1)C=C3/C(=O)NC4=C3C=C(F)C=C4">(mol); // FIXME: stereo
    validate_contains<"CC1=C(C(C2CC2)C3=CC(NS([O-])([O-])C4=CC=C(C=C4)[N+]([O-])=O)=CC=C3)C(=O)OC5=C1CCCCCC5">(mol);
    validate_contains<"CC1=C(C)C(C(=O)C2=CC=C(C)C=C2)=C(N)S1">(mol);
    validate_contains<"CC1=C(C)C(C(=O)C2=CC=CC=C2)=C(N)S1">(mol);
    validate_contains<"CC1=C(C)C(C)=C(C)N1">(mol);
    validate_contains<"CC1=C(C)C(N)=NC=C1">(mol);
    validate_contains<"CC1=C(C)C(O)=CC=C1O">(mol);
    validate_contains<"CC1=C(C)C(O)=NC=C1">(mol);
    validate_contains<"CC1=C(C)C2=C(N1)C=C(C)C(C)=C2">(mol);
    validate_contains<"CC1=C(C)C=C(*)C(=O)N1">(mol);
    validate_contains<"CC1=C(C)C=C(C=C1)C2CCNCC2">(mol);
    validate_contains<"CC1=C(C)C=CN1">(mol);
    validate_contains<"CC1=C(C)N2C=CN=CC2=N1">(mol);
    validate_contains<"CC1=C(C=CC=C1)N2C(=O)C3=C(C)C=CC=C3N=C2CN4C=NC5=C4N=CN=C5N">(mol);
    validate_contains<"CC1=C(CCOP([O-])(=O)OP(O)([O-])=O)SC=[N+]1CC2=CN=C(C)N=C2N">(mol);
    validate_contains<"CC1=C(CCl)N=C(Cl)C(=C1)C#N">(mol);
    validate_contains<"CC1=C(F)C(Cl)=NC(Cl)=N1">(mol);
    validate_contains<"CC1=C(N)C2=C(C=C(F)C(=C2)N3CCOCC3)N=C1">(mol);
    validate_contains<"CC1=C(N)C=CC(N)=C1">(mol);
    validate_contains<"CC1=C(N)C=CC(NC=O)=C1">(mol);
    validate_contains<"CC1=C(N)C=CC(O)=C1">(mol);
    validate_contains<"CC1=C(N)C=CC2=C1N=C(N)C3=CN=CC=C23">(mol);
    validate_contains<"CC1=C(N)N2C=CC=CC2=N1">(mol);
    validate_contains<"CC1=C(N)OC(=N1)C2=CC=CC=C2">(mol);
    validate_contains<"CC1=C(N2C(SC1)C(NC(=O)Cc3ccccc3)C2=O)C(O)=O">(mol);
    validate_contains<"CC1=C(NC2=NC(=CC=N2)C3=CC=CN=C3)C=C(C=C1)C(=O)NC4=CC=CC=C4">(mol);
    //validate_contains<"CC1=C(NC=C1)C=C2/C(=O)NC3=CC=CC=C23">(mol); // FIXME: stereo
}
