#include "Validate.hpp"

void BindingDB_substructure_part_5(OpenBabel::OBMol &mol)
{
    // SMARTS 201 - 250
    validate_contains<"CN(C)CCOC(C1=CC=CC=C1)C1=CC=CC=C1">(mol);
    //validate_contains<"CN(C)C[C@H](O)COc1ccc(Nc2ncc(Br)c(Nc3ccccc3)n2)cc1">(mol); // FIXME: stereo
    validate_contains<"CN(C)S(=O)(=O)C1=CC2=C(C=C1)C=C(Cl)C=C2">(mol);
    validate_contains<"CN(C)S(=O)(=O)C1=CC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"CN(C1=CC2=NN(C)C(C)=C2C=C1)C1=CC=NC(NC2=CC(=C(C)C=C2)S(N)(=O)=O)=N1">(mol);
    validate_contains<"CN(CC1=NC2=C(N=C1)N=C(N)N=C2N)C1=CC=C(C=C1)C(=O)NC(CCC(O)=O)C(O)=O">(mol);
    validate_contains<"CN(CCC(O)C(N1CCC=C1)C1=CCOC=C1)C(O)=O">(mol);
    validate_contains<"CN1C(=O)C(=C)C2=C(C=CC=C2)C1=O">(mol);
    validate_contains<"CN1C(=O)C(=CC2=CN=C(N)N=C12)C1=CC=CC=C1">(mol);
    validate_contains<"CN1C(=O)C(=CC2=CN=C(NC3=CC=CC=C3)N=C12)C1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"CN1C(=O)C(=CC2=CN=CN=C12)C1=CC=CC=C1">(mol);
    validate_contains<"CN1C(=O)C(C)=C(C)C2=CN=C(N)N=C12">(mol);
    validate_contains<"CN1C(=O)C(C)=CC2=CN=C(C)N=C12">(mol);
    validate_contains<"CN1C(=O)C2=CC=CC=C2N=C1C=C">(mol);
    validate_contains<"CN1C(=O)C=CC2=CN=CN=C12">(mol);
    validate_contains<"CN1C(=O)CC2=C(C=CC=C2)C1=O">(mol);
    validate_contains<"CN1C(NC2=CC=C(C=C2)C(F)(F)F)=NC2=CC(OC3=CC(=NC=C3)C3=NC=C(N3)C(F)(F)F)=CC=C12">(mol);
    validate_contains<"CN1C2=CC=CC=C2C2=C1C(=O)NC=C2C1=CC=CC=C1">(mol);
    validate_contains<"CN1C=C(C(C(C)=O)=C1C)C1=CC=CC=C1">(mol);
    validate_contains<"CN1C=C(C)C(C)=C1C">(mol);
    validate_contains<"CN1C=CC(=O)NC1=O">(mol);
    validate_contains<"CN1C=NC2=C(N)N=C(N)N=C12">(mol);
    validate_contains<"CN1CCC(N)(CC1)C(O)=O">(mol);
    validate_contains<"CN1CCN(CC1)C(C1=CC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC2=C(NC(NC3=C(C=CC=C3)C(F)(F)F)=N2)C=C1C(=O)NC1=CC2=C(C=NN2)C=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC2=C(NC(NC3=CC=CC=C3)=N2)C=C1C(=O)NC1=CC2=C(C=NN2)C=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC2=C(NC=N2)C=C1C(=O)NC1=CC2=C(C=NN2)C=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC=C(NC2=NC=C(C)C(NC3=CC(=CC=C3)S(=O)(=O)NC(C)(C)C)=N2)C=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC=C(NC2=NC=C(C)C(NC3=CC=CC=C3)=N2)C=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC=C(NC2=NC=C3C=C(C(=O)N(C)C3=N2)C2=C(Cl)C=CC=C2Cl)C=C1">(mol);
    validate_contains<"CN1CCN(CC1)C1=CC=C2[N]C(=NC2=C1)C1=C2C=CC=CC=C2NC1=O">(mol);
    validate_contains<"CN1CCN(CC1)c1cccc2nc([nH]c12)-c1n[nH]c2cc(ccc12)-c1ccc(N)cc1">(mol);
    validate_contains<"CN1CCNCC1">(mol);
    //validate_contains<"CN1N=C(C(O)=O)C(=N/O)C2=C1C=CC(CC1=CC=C(F)C=C1)=N2">(mol); // FIXME: stereo
    //validate_contains<"CN1N=C(C)C(=N/O)C2=NC=CC=C12">(mol); // FIXME: stereo
    //validate_contains<"CN1N=C(NC(=O)CC2=CC=C(F)C=C2)C(=N/O)C2=NC=CC=C12">(mol); // FIXME: stereo
    //validate_contains<"CN=C(/NCc1ccc(Cl)nc1)N[N+]([O-])=O">(mol); // FIXME: stereo
    validate_contains<"CN=C=N">(mol);
    validate_contains<"CN=CC1=CC=CC=C1O">(mol);
    validate_contains<"CNC(=N)NC">(mol);
    validate_contains<"CNC(=O)C1=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=CC=N1">(mol);
    validate_contains<"CNC(=O)C1=CC2=C(C=C1NC)N=C(N)N2">(mol);
    validate_contains<"CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1">(mol);
    validate_contains<"CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC=CC=C3)C=C2)=C1">(mol);
    validate_contains<"CNC(=O)OC1=CC=CC(O)=C1O">(mol);
    validate_contains<"CNC(=O)OC1=CC=CC2=C1OC(C)(C)O2">(mol);
    validate_contains<"CNC(=O)OC1=CC=CC=C1">(mol);
    //validate_contains<"CNC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC(C)C)CC(=O)NO">(mol); // FIXME: stereo
    validate_contains<"CNC1=C(C(C)C2=C(C=CC=C2)C1=C)N1CCCCC1">(mol);
    validate_contains<"CNC1=C2C=C(O)C=CC2=NC2=C1C=CC=C2">(mol);
}
