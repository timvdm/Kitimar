#include "Validate.hpp"

void BindingDB_substructure_part_2(OpenBabel::OBMol &mol)
{
    // SMARTS 51 - 100
    validate_contains<"CC(=O)C1=C(O)C=C(O)C=C1">(mol);
    validate_contains<"CC(=O)C1=CC=C(C=C1)N1CCN(CC2=CC=CC=C2)CC1">(mol);
    validate_contains<"CC(=O)C1=CCCCC1=O">(mol);
    validate_contains<"CC(=O)C1=NCCS1">(mol);
    validate_contains<"CC(=O)NC(N)=N">(mol);
    validate_contains<"CC(=O)NC1=CC=NC(N)=C1">(mol);
    validate_contains<"CC(=O)NC1=NC(=CO1)C1=NC=CC=C1">(mol);
    validate_contains<"CC(=O)NC1=NC(=CS1)C1=NC=C(C=C1)C(N)=O">(mol);
    validate_contains<"CC(=O)NC1=NC=C(Cl)S1">(mol);
    validate_contains<"CC(=O)NC1=NC=C(SCC2=NC=CO2)S1">(mol);
    validate_contains<"CC(=O)NC1=NC=CN1">(mol);
    validate_contains<"CC(=O)NC1=NC=CN1S(C)(=O)=O">(mol);
    validate_contains<"CC(=O)NC1=NC=CN1S(N)(=O)=O">(mol);
    validate_contains<"CC(=O)NC1=NC=CS1">(mol);
    validate_contains<"CC(=O)NN=CC1=CC=CC=C1">(mol);
    validate_contains<"CC(C)(C)C1=NN(C(NC(=O)NC2=CC=CC=C2)=C1)C1=CC=CC=C1">(mol);
    validate_contains<"CC(C)(NC(CC(=O)NC1=CC=CC=C1)CC1=CC=CC=C1)P(O)(O)=O">(mol);
    validate_contains<"CC(C)=NO">(mol);
    //validate_contains<"CC(C)C1=CN=C(NC(=O)[C@@H](C)C2=CC=C(NC(C)=O)C=C2)S1">(mol); // FIXME: stereo
    validate_contains<"CC(C)N1C=NC2=C(N)N=C(N)N=C12">(mol);
    validate_contains<"CC(C)N1C=NC2=C(NC3=CC=CC=C3)N=C(N)N=C12">(mol);
    validate_contains<"CC(CCCN)NC1=C2N=C(O)C=C(C)C2=C(OC2=CC(=CC=C2)C(F)(F)F)C(O)=C1">(mol);
    validate_contains<"CC(CN)(CS)N1CNC=N1">(mol);
    validate_contains<"CC(N)(CS)N1CNC=N1">(mol);
    validate_contains<"CC(N)(CS)N1N=CNC1=O">(mol);
    validate_contains<"CC(N)=O">(mol);
    validate_contains<"CC(N1CCNCC1)C1=NC2=C(C=CC=C2)C(OC2CCCCC2)=N1">(mol);
    //validate_contains<"CC.OC(CSc1ccccc1)CN1CCC(CC1)Nc1nccc(n1)C(F)(F)F">(mol); // FIXME: components
    validate_contains<"CC1(C)CC2=C3NC(N)=NC3=CC(C(N)=O)=C2O1">(mol);
    validate_contains<"CC1(C)CC2=C3NC(N)=NC3=CC=C2O1">(mol);
    //validate_contains<"CC12CCC3C(C1CCC2=O)[C@H](CC=C)CC1=CC(=O)C=CC31C">(mol); // FIXME: stereo
    validate_contains<"CC12CCC3C(CCC4=C(O)C(=O)CCC34C)C1CCC2=O">(mol);
    validate_contains<"CC1=*OC(=N1)C1=CNC2=CC=CC=C12">(mol);
    validate_contains<"CC1=C(C(=O)N2CCN(CC2)C2=C(Cl)C=C(C=C2)[N+]([O-])=O)C(=NO1)C1=CC=CC=C1">(mol);
    //validate_contains<"CC1=C(C)C(C=C2/C(=O)NC(=O)C3=C2C=CC=C3)=CC=N1">(mol); // FIXME: stereo
    validate_contains<"CC1=C(C)C2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"CC1=C(C)C2=C(O1)C(C)=NC=N2">(mol);
    validate_contains<"CC1=C(C)C2=CC3=C(OC4=C3CCCC4)C=C2OC1=O">(mol);
    validate_contains<"CC1=C(C=NN1)C1=CC=CC=C1">(mol);
    //validate_contains<"CC1=C(CCC(=O)N[C@H](CSCC2=CC=CC=C2)C(=O)NCCCC(O)=O)C(=O)OC2=CC3=C(C=C12)C1=C(CCCC1)O3">(mol); // FIXME: stereo
    validate_contains<"CC1=C(Cl)C=CC=C1">(mol);
    validate_contains<"CC1=C(NC2=CC=CC=C2C(O)=O)C=CC=C1Cl">(mol);
    validate_contains<"CC1=C(O)C(=O)C2=C(O)C=C(O)C(C)=C2O1">(mol);
    validate_contains<"CC1=C(O)C(=O)C2=C(O)C=C(O)C=C2O1">(mol);
    validate_contains<"CC1=C(O)C(=O)C2=CC=C(O)C(C)=C2O1">(mol);
    validate_contains<"CC1=CC(=NC2=NC=NC(N)=C12)C1=CC=C(N)C=C1">(mol);
    validate_contains<"CC1=CC(=NC2=NC=NC(N)=C12)C1=CN=C(N)C=C1">(mol);
    validate_contains<"CC1=CC(=O)C2=C(O)C=C(O)C(C)=C2O1">(mol);
    validate_contains<"CC1=CC(=O)C2=CC=C(O)C(C)=C2O1">(mol);
    validate_contains<"CC1=CC(=O)NN1">(mol);
}
