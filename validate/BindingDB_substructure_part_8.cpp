#include "Validate.hpp"

void BindingDB_substructure_part_8(OpenBabel::OBMol &mol)
{
    // SMARTS 351 - 400
    validate_contains<"NC1=C(SC=N1)C(=O)C1=CC=CC=C1">(mol);
    validate_contains<"NC1=C2C(=CC(=NC2=NC=N1)C1=CN=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"NC1=C2C=C(O)C=CC2=NC2=C1C=CC=C2">(mol);
    validate_contains<"NC1=C2C=CC=CC2=NC2=C1C=CC=C2">(mol);
    validate_contains<"NC1=C2N=C(O)C=CC2=C(OC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"NC1=C2N=CC=CC2=C(OC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"NC1=CC(=NC2=CC=C(C=C2)S(=O)=O)C2=C(C=CC=C2)C1=O">(mol);
    validate_contains<"NC1=CC(=NC2=CC=CC=C2)C2=C(C=CC=C2)C1=O">(mol);
    validate_contains<"NC1=CC(NC2=NC=CC(=N2)C2=CC=CN=C2)=CC=C1">(mol);
    validate_contains<"NC1=CC2=C(C=C1)N=CC(C#N)=C2N">(mol);
    validate_contains<"NC1=CC2=C(SC3=C2N=CN=C3NC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"NC1=CC=C(C=C1)C1=NC(=C(N1)C1=CC=NC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"NC1=CC=C(C=N1)C1=NC(NC2=CC=CC=C2)=NC=C1">(mol);
    //validate_contains<"NC1=CC=C2NC(=O)C(C2=C1)=C(/NC1=CC=C(CN2CCCCC2)C=C1)C1=CC=CC=C1">(mol); // FIXME: stereo
    validate_contains<"NC1=CC=CC2=C1C(S)=CC=C2">(mol);
    validate_contains<"NC1=CC=CC2=C1N=CC=C2">(mol);
    validate_contains<"NC1=NC(=CC(C2=CC=CC=C2)=C1C#N)C1=C(O)C=CC=C1">(mol);
    validate_contains<"NC1=NC(=CC=N1)C1=CC=CC2=C1NCCN2">(mol);
    validate_contains<"NC1=NC(=CN1)C1=CC=NC=C1">(mol);
    validate_contains<"NC1=NC(=CN1)C1=NC=CC=C1">(mol);
    validate_contains<"NC1=NC(=CO1)C1=CC=CC=C1">(mol);
    validate_contains<"NC1=NC(=CS1)C1=NC(=CC=C1)S(N)(=O)=O">(mol);
    validate_contains<"NC1=NC(=CS1)C1=NC=CC=C1">(mol);
    validate_contains<"NC1=NC(=CS1)N1C=CC=N1">(mol);
    validate_contains<"NC1=NC(=O)N(C=C1)C2OC(COP(O)(O)=O)C(O)C2O">(mol);
    validate_contains<"NC1=NC(NC2=CC=CC=C2)=NN1C(=O)C1=CC=CC=C1">(mol);
    validate_contains<"NC1=NC(NCC2=CC=CC=C2)=C2N=CNC2N1">(mol);
    validate_contains<"NC1=NC2=C(N=CN2COCCO)C(=O)N1">(mol);
    validate_contains<"NC1=NC2=C(S1)C(CC(=O)N2)C1=CC(Cl)=CC=C1">(mol);
    validate_contains<"NC1=NC2=CC=C3OCCC3=C2N1">(mol);
    validate_contains<"NC1=NC=C(C2=CC=CC=C2)C2=C1C=C(C=C2)C(O)=O">(mol);
    validate_contains<"NC1=NC=C(S1)C1=NC(NC2=CC=CC=C2)=NC=C1">(mol);
    validate_contains<"NC1=NC=C2C=C(C(=O)N(O)C2=N1)C1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"NC1=NC=C2C=C(C(=O)NC2=N1)C1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"NC1=NC=C2N=CNC2N1">(mol);
    validate_contains<"NC1=NC=CC(=C1)C1=NC(NC2=CC=CC(Cl)=C2)=NC=C1">(mol);
    validate_contains<"NC1=NC=CC2=C1C=C(C=C2)C(O)=O">(mol);
    validate_contains<"NC1=NC=NC2=C1C=C(C=C2)C(O)=O">(mol);
    validate_contains<"NC1=NC=NC2=C1N=CN2C1CCCO1">(mol);
    validate_contains<"NC1=NNC2=CC=C(N)C=C12">(mol);
    validate_contains<"NC1C=C(N)C(=O)C2=C1C=CC=C2">(mol);
    validate_contains<"NC1CC2N(CCC3=CC=CC=C23)CC1N1CCCC1=O">(mol);
    validate_contains<"NC1CCCC1C(=O)N1CCCC1">(mol);
    validate_contains<"NC1CCCC1C1=CC=CC=C1">(mol);
    validate_contains<"NC1CCCCC1C1=CC=CC=C1">(mol);
    validate_contains<"NC1CCCCC1N1CCCC1=O">(mol);
    validate_contains<"NC1CCNCC1N1CCCC1=O">(mol);
    validate_contains<"NC1CNNC1=O">(mol);
    //validate_contains<"NC=C1/C(=O)NC(=O)C2=C1C=CC=C2">(mol); // FIXME: stereo
    validate_contains<"NCC(=O)NC1=CC(NC2=NC=CC(=N2)C2=CC=CN=C2)=CC=C1">(mol);
}
