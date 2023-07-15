#include "Validate.hpp"

void BindingDB_substructure_part_3(OpenBabel::OBMol &mol)
{
    // SMARTS 101 - 150
    //validate_contains<"CC1=CC(C)=C(N1)C=C1/C(=O)NC2=CC=C(F)C=C12">(mol); // FIXME: stereo
    validate_contains<"CC1=CC(C)=NC(=N1)S(N)(=O)=O">(mol);
    validate_contains<"CC1=CC(C2NC3=C(N2)C=CC=C3)=C2C=CC=CC(O)=C12">(mol);
    //validate_contains<"CC1=CC(Cl)=C(N1)C=C1/C(=O)NC2=CC=C(F)C=C12">(mol); // FIXME: stereo
    validate_contains<"CC1=CC(NC2=CC=CC=C2)=NC=N1">(mol);
    validate_contains<"CC1=CC(NC2=NC=NC(C)=C2)=CC=C1">(mol);
    validate_contains<"CC1=CC(O)=NC2=C(N)C=C(O)C(OC3=CC(=CC=C3)C(F)(F)F)=C12">(mol);
    validate_contains<"CC1=CC(O)=NC2=C(N)C=C(O)C(OC3=CC=CC=C3)=C12">(mol);
    validate_contains<"CC1=CC2=NC=CN2C(N)=N1">(mol);
    validate_contains<"CC1=CC=C(C=C1)C1=NC(=C(N1)C1=CC=NC=C1)C1=CC=C(F)C=C1">(mol);
    validate_contains<"CC1=CC=C(C=C1)N1N=C(C=C1NC(=O)NC1=CC=C(OCCC2=CC=CN=C2)C2=CC=CC=C12)C(C)(C)C">(mol);
    validate_contains<"CC1=CC=C(C=C1)N1N=C(C=C1NC(=O)NC1=CC=C(OCCC2=CC=NC=C2)C2=CC=CC=C12)C(C)(C)C">(mol);
    validate_contains<"CC1=CC=C(C=C1)N1N=C(C=C1NC(=O)NC1=CC=C(OCCN2C=CN=C2)C2=CC=CC=C12)C(C)(C)C">(mol);
    validate_contains<"CC1=CC=C(C=C1)S(=O)(=O)N([*])C([*])C(=O)NO">(mol);
    //validate_contains<"CC1=CC=C(N1)C=C1/C(=O)NC2=CC=C(F)C=C12">(mol); // FIXME: stereo
    validate_contains<"CC1=CC=C(O)C(=C1)C1=C(NN=C1)C1=CC=CC=C1">(mol);
    validate_contains<"CC1=CC=C(O)C2=C1C(=CN=C2)C1CCC(O)C1">(mol);
    validate_contains<"CC1=CC=C(O)C2=C1N1C=CC(O)=C1CC2">(mol);
    validate_contains<"CC1=CC=CC(=N1)S(N)(=O)=O">(mol);
    validate_contains<"CC1=CC=CC(C2=CC=CC=C2)=C1CN">(mol);
    validate_contains<"CC1=CC=CC(S)=N1">(mol);
    validate_contains<"CC1=CC=NC(=N1)S(N)(=O)=O">(mol);
    validate_contains<"CC1=CC=NC(S)=N1">(mol);
    validate_contains<"CC1=CN2C(N=CC=C2C)=C1C">(mol);
    validate_contains<"CC1=CN=C(C=C1C)C1=CSC(N)=N1">(mol);
    validate_contains<"CC1=CN=C(N)C2=C1C=CC(=C2)C(O)=O">(mol);
    validate_contains<"CC1=CN=C(N=C1)S(N)(=O)=O">(mol);
    validate_contains<"CC1=CN=C(S)N=C1">(mol);
    validate_contains<"CC1=COC2=C1C(=O)C(=O)C1=C2C=CC2=C1C=CC=C2C">(mol);
    validate_contains<"CC1=NC2=C(C(=O)N1)C(F)=CC=C2">(mol);
    validate_contains<"CC1=NC2=C(C=CC=C2)N1C(=O)C1=CC(Cl)=CC=C1Cl">(mol);
    validate_contains<"CC1=NC2=C(C=CC=C2)N1C(=O)C1=CC=CC=C1">(mol);
    validate_contains<"CC1=NC2=NC=NC(N)=C2C(C)=C1">(mol);
    validate_contains<"CC1=NC=C2C=C(C)C(=O)N(C3CCCC3)C2=N1">(mol);
    validate_contains<"CC1=NC=NC(C)=N1">(mol);
    validate_contains<"CC1=NN=C(O1)C1=C(C=CS1)S(=O)=O">(mol);
    validate_contains<"CC1=NN=C(O1)C1=C(C=CS1)S(C)(=O)=O">(mol);
    validate_contains<"CC1=NN=C(O1)C1=C(C=CS1)S(N)(=O)=O">(mol);
    validate_contains<"CC1=NN=C(O1)C1=C(S)C=CS1">(mol);
    validate_contains<"CC1=NN=C2N1C1=C(C=CC=C1)N=C2C1=CC(Cl)=NN1">(mol);
    validate_contains<"CC1=NNC(=C1)C1=NC2=C(C=CC=C2)N2C(C)=NN=C12">(mol);
    validate_contains<"CC1CC(C)C(N)C1">(mol);
    validate_contains<"CC1CC(C)NN1">(mol);
    validate_contains<"CC1CC=C1">(mol);
    validate_contains<"CC1CCC(C)C1N">(mol);
    validate_contains<"CC1CCC(C)C1NC(N)N">(mol);
    validate_contains<"CC1CCC(C)N1C(N)=N">(mol);
    validate_contains<"CC1CCC(C1N)C1=CC=CC=C1">(mol);
    validate_contains<"CC1CCCC2=CC=CC=C12">(mol);
    validate_contains<"CC1CCCCC1">(mol);
}
