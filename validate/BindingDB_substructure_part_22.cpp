#include "Validate.hpp"

void BindingDB_substructure_part_22(OpenBabel::OBMol &mol)
{
    // SMARTS 1051 - 1100
    validate_contains<"CC1=C(NN=C1)C2=CNC=C2">(mol);
    validate_contains<"CC1=C(O)C(=CC=C1)C2=CC3=C(N2)C=CC(=C3)C(N)=N">(mol);
    validate_contains<"CC1=C(O)C=C2C=CC(N)=CC2=C1">(mol);
    //validate_contains<"CC1=C(O)C=CC(C=C2/CCCC2=O)=C1">(mol); // FIXME: stereo
    validate_contains<"CC1=C(O)C=CC(C=CC(=O)CC(=O)C=CC2=CC(C)=C(O)C=C2)=C1">(mol);
    validate_contains<"CC1=C(P)[N+](C)=CC=C1">(mol);
    validate_contains<"CC1=C(SC(=N1)C2COC3=C(O2)C=CC=C3)C(O)=O">(mol);
    //validate_contains<"CC1=C(SC(=N1)[C@@H]2COC3=C(O2)C=CC=C3)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"CC1=C(SC(N)=N1)C2=CC(NC3=CC=CC(=C3)N(=O)=O)=CC=C2">(mol);
    validate_contains<"CC1=C(SC(N)=N1)C2=NC(NC3=CC(=CC=C3)N(=O)=O)=NC=C2">(mol);
    validate_contains<"CC1=C(SC=N1)C2=CC=NC(N)=N2">(mol);
    validate_contains<"CC1=C2N=CC3=C(NC(=O)NC3=O)C2=CC=C1">(mol);
    validate_contains<"CC1=C2N=CNC2=NC(N)=N1">(mol);
    validate_contains<"CC1=CC(=CC(=C1)C(O)=O)C(O)=O">(mol);
    validate_contains<"CC1=CC(=CC(C)=C1)C(=O)C2=C(C)C(=O)NC(C)=C2C">(mol);
    validate_contains<"CC1=CC(=CC(C)=C1)C(=O)C2=C(C)C(C)=NC(=O)C2N">(mol);
    validate_contains<"CC1=CC(=CC=C1N)C2=NC3=C(S2)C=CC=C3">(mol);
    validate_contains<"CC1=CC(=CS1)C2=CC=CC=C2">(mol);
    validate_contains<"CC1=CC(=O)C(O)=C(CN2CCCCC2)O1">(mol);
    validate_contains<"CC1=CC(=O)C2=C(O)C=CC=C2C1=O">(mol);
    validate_contains<"CC1=CC(=O)OC2=C1C=C(OCC(O)=O)C(OCC(O)=O)=C2">(mol);
    validate_contains<"CC1=CC(=O)OC2=C1C=CC(NC(=O)C3=CC=CC=C3)=C2">(mol);
    validate_contains<"CC1=CC(Br)=CN=C1N">(mol);
    validate_contains<"CC1=CC(C(N)=O)=C(C)N1">(mol);
    validate_contains<"CC1=CC(C)(C)NC(=O)N1">(mol);
    validate_contains<"CC1=CC(C)=C(C)O1">(mol);
    validate_contains<"CC1=CC(C)=CC(C)=C1">(mol);
    validate_contains<"CC1=CC(C)=CC=C1">(mol);
    validate_contains<"CC1=CC(C)=NC(N)=N1">(mol);
    validate_contains<"CC1=CC(C)=NC=N1">(mol);
    validate_contains<"CC1=CC(CS1)C2=CC=CC=C2">(mol);
    validate_contains<"CC1=CC(N)=CC2=C1NC(=N2)C3=C(N)C=CNC3=O">(mol);
    validate_contains<"CC1=CC(NC(=O)C2=CC=CC=C2)=NN1">(mol);
    validate_contains<"CC1=CC(NC2=CC=NC(NC3=CC=C(N)C=C3)=N2)=NN1">(mol);
    validate_contains<"CC1=CC(NC2=NC(=CC=N2)N3C=CN=C3C4=CC=CC=C4)=CC(C)=C1">(mol);
    validate_contains<"CC1=CC(NC2=NC(=CC=N2)N3C=CN=C3C4=CN=CC=C4)=CC(C)=C1">(mol);
    validate_contains<"CC1=CC(O)=CC(Cl)=C1">(mol);
    validate_contains<"CC1=CC(OC(C)(C)C)=CN1">(mol);
    validate_contains<"CC1=CC(S)=NC=C1">(mol);
    validate_contains<"CC1=CC2=C(C=C1)C(Cl)=CC(=C2)C3=CC(Cl)=CC(C)=C3">(mol);
    validate_contains<"CC1=CC2=C(C=C1)N=CC=C2">(mol);
    validate_contains<"CC1=CC2=C(C=CC=C2)C(C)=N1">(mol);
    validate_contains<"CC1=CC2=C(C=CC=C2)C(N)=N1">(mol);
    validate_contains<"CC1=CC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"CC1=CC2=C(C=CC=C2)C=C1C">(mol);
    validate_contains<"CC1=CC2=C(C=NN2)C=C1">(mol);
    validate_contains<"CC1=CC2=C(N1)N=CC(C)=C2">(mol);
    validate_contains<"CC1=CC2=C(NC=N2)C=C1">(mol);
    validate_contains<"CC1=CC2=C(NN=C2)C=C1">(mol);
    validate_contains<"CC1=CC2=C(NN=C2C)C=C1">(mol);
}
