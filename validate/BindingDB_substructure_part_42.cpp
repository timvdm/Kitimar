#include "Validate.hpp"

void BindingDB_substructure_part_42(OpenBabel::OBMol &mol)
{
    // SMARTS 2051 - 2100
    validate_contains<"NC1=CC2=C(NC=N2)C=C1">(mol);
    validate_contains<"NC1=CC2=C(OCO2)C=C1">(mol);
    validate_contains<"NC1=CC=C(C=C1)C(=O)C2=CC=CC=C2">(mol);
    validate_contains<"NC1=CC=C(C=C1)C1=NC(=C(N1)C1=CC=NC=C1)C1=CC=C(F)C=C1">(mol);
    validate_contains<"NC1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"NC1=CC=C(C=C1)S(=O)(=O)NC2=CC=CC=C2">(mol);
    //validate_contains<"NC1=CC=C(C[C@@H]2[C@H](O)[C@@H](O)[C@@H](CC3=CC=C(N)C=C3)N(CC4=CC=CC=C4)C(=O)N2CC5=CC=CC=C5)C=C1">(mol); // FIXME: stereo
    validate_contains<"NC1=CC=C(NC(O)CCC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"NC1=CC=C2C(=C1)C=C(C(N=NC3=CC=C4C(O)=CC(=CC4=C3)S(O)(O)O)=C2O)S(O)(O)O">(mol);
    validate_contains<"NC1=CC=C2C=C(N)C=NC2=C1">(mol);
    validate_contains<"NC1=CC=C2NC3=C4C=CC=CC4=NC=C3C(N)=NC2=C1">(mol);
    validate_contains<"NC1=CC=CC(=C1)C2=CC=CC=C2">(mol);
    validate_contains<"NC1=CC=CC(Cl)=C1">(mol);
    validate_contains<"NC1=CC=CC2=C1C=CC=C2[S](=O)=O">(mol);
    validate_contains<"NC1=CC=CC2=C1C=CC=C2[S](O)O">(mol);
    validate_contains<"NC1=CC=CN=C1S">(mol);
    validate_contains<"NC1=CC=NC2=C1C=C(O)C(O)=C2">(mol);
    validate_contains<"NC1=CC=NC2=C1N=CN2C3OC(O)CC3O">(mol);
    validate_contains<"NC1=CNC(=O)C2=C1SC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"NC1=NC(=CC(C2=CC(NC=O)=CC=C2)=C1C#N)C1=C(O)C=CC=C1">(mol);
    validate_contains<"NC1=NC(=CO1)C1=NC=CC=C1">(mol);
    validate_contains<"NC1=NC(=CS1)C2CCNCC2">(mol);
    validate_contains<"NC1=NC(=O)CC=N1">(mol);
    validate_contains<"NC1=NC(N)=C(CC2=CC=CC=C2)C=N1">(mol);
    validate_contains<"NC1=NC(N)=C(S1)C(=O)C1=CC=CC=C1">(mol);
    validate_contains<"NC1=NC(N)=C2C=C(CCC3=CC=C(C=C3)C(=O)NCCC(O)CC(O)CC(O)=O)C=CC2=N1">(mol);
    validate_contains<"NC1=NC(N)=C2N=C(C3=CC=CC(O)=C3)C(=NC2=N1)C4=CC(O)=CC=C4">(mol);
    validate_contains<"NC1=NC(N)=NC(Cl)=C1">(mol);
    validate_contains<"NC1=NC(S)=NC(Cl)=C1">(mol);
    validate_contains<"NC1=NC(S)=NC=C1">(mol);
    validate_contains<"NC1=NC2=C(C(=O)N1)C(CCNNCC3=CN(C4OC(CO)CC4O)C5=C3C(=O)NC(N)=N5)=CN2C6CC(O)C(CO)O6">(mol);
    validate_contains<"NC1=NC2=C(C=CC(=O)N2)C=N1">(mol);
    validate_contains<"NC1=NC2=C(C=CC=C2)C3=NC=NN13">(mol);
    validate_contains<"NC1=NC2=C(C=CC=C2N1)C(=O)NC3=NCCN3">(mol);
    validate_contains<"NC1=NC2=C(C=CN2)C=N1">(mol);
    validate_contains<"NC1=NC2=C(N1)C=C(C=O)C(N)=C2">(mol);
    validate_contains<"NC1=NC2=C(N1)C=CC(N)=N2">(mol);
    validate_contains<"NC1=NC2=C(N=CN2)C(=O)N1">(mol);
    validate_contains<"NC1=NC2=C(N=NN2)C(=O)N1">(mol);
    validate_contains<"NC1=NC2=C(S1)C(N)=NC=N2">(mol);
    validate_contains<"NC1=NC2=CC=C(N)N=C2S1">(mol);
    validate_contains<"NC1=NC2=CN=CN=C2S1">(mol);
    validate_contains<"NC1=NC=C(Br)C=N1">(mol);
    validate_contains<"NC1=NC=C(S)S1">(mol);
    validate_contains<"NC1=NC=C(S1)N(=O)=O">(mol);
    validate_contains<"NC1=NC=C2C=CC=CC2=N1">(mol);
    validate_contains<"NC1=NC=CC(=C1)C1=NC(NC2=CC=CC=C2)=NC=C1">(mol);
    validate_contains<"NC1=NC=CC(=N1)N2C=NC3=CC=CC=C23">(mol);
    validate_contains<"NC1=NC=CC(NC2=NOC=C2)=N1">(mol);
    validate_contains<"NC1=NC=NC2=C1C(C3=CC=CC=C3)=NN2C">(mol);
}
