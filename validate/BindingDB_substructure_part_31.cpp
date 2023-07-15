#include "Validate.hpp"

void BindingDB_substructure_part_31(OpenBabel::OBMol &mol)
{
    // SMARTS 1501 - 1550
    validate_contains<"CN1C(=O)N(C2CCCC2)C3=CC=CC=C3C1=O">(mol);
    validate_contains<"CN1C(=O)N2CCC3=C(NC4=C3C=CC=C4)C2(C)C1=O">(mol);
    validate_contains<"CN1C(=O)NC2=CC=CC=C2C1=O">(mol);
    //validate_contains<"CN1C(=S)[C@@H](C(=O)NC2=CC=CC=C2)C3=CC=CC=C13">(mol); // FIXME: stereo
    validate_contains<"CN1C(C)=C(C)C(=O)N(C)C1=O">(mol);
    validate_contains<"CN1C(C)=C(C)C2=C1C=CC(C)=C2">(mol);
    validate_contains<"CN1C(C)=NN=C1C">(mol);
    validate_contains<"CN1C(CCCC1CC(O)C2=CC3=C(C=CC=C3)C=C2)CC(O)C4=CC5=CC=CC=C5C=C4">(mol);
    validate_contains<"CN1C(CCCC1CC(O)C2=CC=CC=C2)CC(O)C3=CC=CC=C3">(mol);
    validate_contains<"CN1C2=C(C=CC=C2)C(=O)C3=C1C=C(O)C=C3O">(mol);
    validate_contains<"CN1C2=C(N3C=CNC3=N2)C(=O)NC1=O">(mol);
    validate_contains<"CN1C=C(C(C)=C1C)C1=CC=CC=C1">(mol);
    validate_contains<"CN1C=C(C)C2=C1C=CC(C)=C2">(mol);
    validate_contains<"CN1C=C(C=S1)C2=CC=CC=C2">(mol);
    validate_contains<"CN1C=CC2=C1C3=C(C=CN3C)C4=C2CNC4=C">(mol);
    validate_contains<"CN1C=CC2=C3CNC(=O)C3=C4C=CN(C)C4=C12">(mol);
    validate_contains<"CN1C=CC=C1">(mol);
    validate_contains<"CN1C=CC=CC1=O">(mol);
    validate_contains<"CN1C=CN=CC1=O">(mol);
    validate_contains<"CN1C=NC(=C1)C2=CC3=NC=CC(OC4=C(F)C=C(NC(=S)NC(=O)CC5=CC=CC=C5)C=C4)=C3S2">(mol);
    validate_contains<"CN1C=NC(=N1)C(N)=O">(mol);
    validate_contains<"CN1C=NC2=C1C(=O)NC(=O)N2C">(mol);
    validate_contains<"CN1C=NC2=C1C3=C(C=CC=C3)N=C2N">(mol);
    validate_contains<"CN1C=NC2=C1C=CC=C2">(mol);
    validate_contains<"CN1C=NC2=CN=CN=C12">(mol);
    validate_contains<"CN1C=NN=C1C">(mol);
    validate_contains<"CN1CC(=O)N(C)CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
    validate_contains<"CN1CC(=O)NC2=CC=CC=C12">(mol);
    validate_contains<"CN1CC(C(N)=O)C2=C1C=CC=C2">(mol);
    validate_contains<"CN1CC(C=C2C1CC3=CNC4=CC=CC2=C34)C(N)=O">(mol);
    validate_contains<"CN1CC(F)(F)CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
    validate_contains<"CN1CC(O)C(O)CN(C)C1=O">(mol);
    validate_contains<"CN1CC2=CC=CC=C2NCC1=O">(mol);
    validate_contains<"CN1CC2=CN=C3C=CC=CC3=C2N(C)C1=O">(mol);
    validate_contains<"CN1CCC(C(O)C1)C(=O)NCC2=CN(C)C=N2">(mol);
    validate_contains<"CN1CCC(CC1)C2=CC(C)=C(C)C=N2">(mol);
    validate_contains<"CN1CCC(CC1)C2=CC=CC=C2">(mol);
    validate_contains<"CN1CCC(CC1)C2=CC=CN=N2">(mol);
    validate_contains<"CN1CCC(CC1)C2N=C(C)CN2C3=CC=CC=C3">(mol);
    validate_contains<"CN1CCC(CC1)OC(C2=CC=C(F)C=C2)C3=CC=C(F)C=C3">(mol);
    validate_contains<"CN1CCC(CC2=CC=CC=C2)C(O)C1">(mol);
    validate_contains<"CN1CCC(N[S](=O)=O)C1=O">(mol);
    validate_contains<"CN1CCC1=O">(mol);
    validate_contains<"CN1CCC2=C(C1)C3=C(C=CC(C)=C3)N2CCC4=CN=C(C)C=C4">(mol);
    validate_contains<"CN1CCC=C(C)C1">(mol);
    validate_contains<"CN1CCCC(O)C1">(mol);
    validate_contains<"CN1CCCN(C2CCCN3C(=O)C(O)=C(N=C23)C(=O)NCc4ccc(F)cc4)S1(=O)=O">(mol);
    validate_contains<"CN1CCN(CC(=O)NC2=CC(NC3=NC=CC(=N3)C3=CC=CN=C3)=CC=C2)CC1">(mol);
    //validate_contains<"CN1CCN(CC1(C)C)[C@@H]2C[C@H](C3=C2C=C(Cl)C=C3)C4=CC=CC=C4">(mol); // FIXME: stereo
    validate_contains<"CN1CCN(CC1)C(=O)C(=O)NC(C)(C)C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
}
