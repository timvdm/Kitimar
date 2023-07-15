#include "Validate.hpp"

void BindingDB_substructure_part_23(OpenBabel::OBMol &mol)
{
    // SMARTS 1101 - 1150
    validate_contains<"CC1=CC2=C(OC(=O)C(=C2)C(=O)OC3=CC=CC=C3)C=C1">(mol);
    validate_contains<"CC1=CC2=C(OCO2)C(N)=C1">(mol);
    validate_contains<"CC1=CC2=C(OCO2)C=C1N">(mol);
    validate_contains<"CC1=CC2=C(OCO2)C=C1N(=O)=O">(mol);
    validate_contains<"CC1=CC2=CC3=C(OCCO3)C=C2CN1">(mol);
    validate_contains<"CC1=CC2=CC=CC=C2CN1">(mol);
    validate_contains<"CC1=CC=C(C)N1C2=CC=CC=C2">(mol);
    validate_contains<"CC1=CC=C(C=C1)C(N)=N">(mol);
    validate_contains<"CC1=CC=C(C=C1)C2=CC(C3=CC=CC=C3)=C(C#N)C(SCC(=O)NC4=CC=C(C=C4)C(O)=O)=N2">(mol);
    validate_contains<"CC1=CC=C(C=C1)N1N=C(C=C1NC(N)=O)C(C)(C)C">(mol);
    validate_contains<"CC1=CC=C(C=C1)N2N=C(C=C2NC(=O)NC3=C4C=CC=CC4=C(OCCN5CCOCC5)C=C3)C(C)(C)C">(mol);
    validate_contains<"CC1=CC=C(C=C1)N2N=CC=C2NC(=O)NC3=C4C=CC=CC4=C(O)C=C3">(mol);
    validate_contains<"CC1=CC=C(C=C1)N=CC=CC2=CC=CC=C2">(mol);
    //validate_contains<"CC1=CC=C(C=C1)S(=O)(=O)NN1C(=S)SC(=C/C2=CC=C(Cl)C=C2Cl)C1=O">(mol); // FIXME: stereo
    validate_contains<"CC1=CC=C(C=C1C)C2CC3CC3C2">(mol);
    validate_contains<"CC1=CC=C(C=C1S(=O)(=O)NC2CCCCC2)C(O)N3CCC4CCCCC4C3">(mol);
    validate_contains<"CC1=CC=C(C=CC(=O)CC(=O)C=CC2=CC(O)=C(C)C=C2)C=C1O">(mol);
    validate_contains<"CC1=CC=C(C=CC(=O)OCCCC2=CC=C(O)C=C2)C=C1">(mol);
    validate_contains<"CC1=CC=C(CCCOC(=O)C=CC2=CC=C(C)C=C2)C=C1">(mol);
    validate_contains<"CC1=CC=C(COC2=C(SC(=C2)N3C=NC4=C3C=CC=C4)C(N)=O)C=C1">(mol);
    //validate_contains<"CC1=CC=C(N1)C=C2/C(=O)NC3=CC=CC=C23">(mol); // FIXME: stereo
    validate_contains<"CC1=CC=C(NS(=O)(=O)C2=CC=CC=C2)C(=C1)C(O)=O">(mol);
    validate_contains<"CC1=CC=C(O1)[N+]([O-])=O">(mol);
    validate_contains<"CC1=CC=C2NC=[N+](C)C2=C1">(mol);
    //validate_contains<"CC1=CC=CC(CN2[C@H](CC3=CC=CC=C3)[C@H](O)[C@@H](O)[C@@H](CC4=CC=CC=C4)N(CC5=CC=CC(C)=C5)C2=O)=C1">(mol); // FIXME: stereo
    validate_contains<"CC1=CC=CC(NC(=O)C2=CC(=CN=C2)C3=CC=CC=C3)=C1">(mol);
    validate_contains<"CC1=CC=CC(NC2=CC(Cl)=NC(SCC(O)=O)=N2)=C1C">(mol);
    validate_contains<"CC1=CC=CC=C1N2C=CC(=C2)N(=O)=O">(mol);
    //validate_contains<"CC1=CC=CN1.CCC2=CNC=C2">(mol); // FIXME: components
    validate_contains<"CC1=CC=CN1C2=CC=CC=C2">(mol);
    validate_contains<"CC1=CC=CS1">(mol);
    validate_contains<"CC1=CC=NC(=C1)C2=NC=C(S2)C3=CC=CC=C3">(mol);
    validate_contains<"CC1=CC=NC(O)=N1">(mol);
    validate_contains<"CC1=CCN(CC2=CC=C(C=C2)C(=O)NC3=CC=NC(NC4=CC=CC(=C4)N5CC=CN=C5)=C3)C=N1">(mol);
    validate_contains<"CC1=CCN2CC(=NC2=C1)C3=CC=C(F)C=C3">(mol);
    validate_contains<"CC1=CN(C2CC2)C3=C(F)C(C)=C(F)C=C3C1=O">(mol);
    validate_contains<"CC1=CN(C=N1)C2=CC(NC(=O)C3=CC=C(C)C(NC4=NC=CC(=N4)C5=CC=CN=C5)=C3)=CC(=C2)C(F)(F)F">(mol);
    validate_contains<"CC1=CN=C(N)C(Br)=C1">(mol);
    validate_contains<"CC1=CN=C(N)C(N)=C1">(mol);
    validate_contains<"CC1=CN=C(N)S1">(mol);
    validate_contains<"CC1=CN=C(S)C=C1">(mol);
    validate_contains<"CC1=CN=C2C(=CC=CC2=C1C3=CC(OCC4=CC=C(CC(O)=O)C=C4)=CC=C3)C(F)(F)F">(mol);
    validate_contains<"CC1=CN=CN=C1">(mol);
    validate_contains<"CC1=CN=NC(N)=N1">(mol);
    validate_contains<"CC1=CNC=C1">(mol);
    //validate_contains<"CC1=CNC=C1.CC2CCCCC2.CC3CCCC4(CCCCC4)C3">(mol); // FIXME: components
    validate_contains<"CC1=NC(=CC(CC2=NC=C(S2)C(=O)NC3=C(Cl)C=CC=C3C)=N1)N4CCNCC4">(mol);
    validate_contains<"CC1=NC(=CC=N1)C2=C(NC(=N2)C3OCC(C)(CO3)C(=O)N4CCOCC4)C5=CC=C(F)C=C5">(mol);
    validate_contains<"CC1=NC(C)=C(N)O1">(mol);
    validate_contains<"CC1=NC(C)=C(O1)N2C=CC=C2">(mol);
}
