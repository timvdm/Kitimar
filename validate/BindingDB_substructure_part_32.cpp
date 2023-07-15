#include "Validate.hpp"

void BindingDB_substructure_part_32(OpenBabel::OBMol &mol)
{
    // SMARTS 1551 - 1600
    validate_contains<"CN1CCN(CC1)C2=CC(Cl)=C(Cl)C=C2">(mol);
    validate_contains<"CN1CCN(CC1)C2=Nc3cc(Cl)ccc3Nc4ccccc24">(mol);
    validate_contains<"CN1CCN(CC1)S(=O)(=O)C2=CC=CC=C2">(mol);
    validate_contains<"CN1CCN(CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C)C(=O)OC(C)(C)C">(mol);
    validate_contains<"CN1CCN(CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C)C(=O)c4ccccc4">(mol);
    validate_contains<"CN1CCN(CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C)C(C)=O">(mol);
    validate_contains<"CN1CCN(CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C)S(C)(=O)=O">(mol);
    validate_contains<"CN1CCN(CC2=CC=C(C=C2)C(=O)NC2=CC(NC3=NC=CC(=N3)C3=CC=CN=C3)=CC=C2)CC1">(mol);
    validate_contains<"CN1CCN(CC2=CC=C(C=C2)C(=O)NC3=CC=C(C)C(NC4=NC=CC(=N4)C5=CC=CN=C5)=C3)CC1">(mol);
    //validate_contains<"CN1CCN(C[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C)S(=O)(=O)c4ccccc4">(mol); // FIXME: stereo
    validate_contains<"CN1CCN(c(cc2)cc3c2nc(c(c(N)c(c(F)ccc4)c4[nH]5)c5=O)[nH]3)CC1">(mol);
    validate_contains<"CN1CCN2C(C1)C3=CC=CC=C3CC4=C2C=CC=C4">(mol);
    validate_contains<"CN1CCN2C(C1)C3=CC=CC=C3CC4=C2N=CC=C4">(mol);
    validate_contains<"CN1CCOCC1C2=NC(C(=O)NCc3ccc(F)c(Br)c3)=C(O)C(=O)N2C">(mol);
    validate_contains<"CN1CCOCC1C2=NC(C(=O)NCc3ccc(F)c(F)c3)=C(O)C(=O)N2C">(mol);
    validate_contains<"CN1CCOCC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2">(mol);
    validate_contains<"CN1CCOCC1C2=NC(C(=O)NCc3cccc(Br)c3)=C(O)C(=O)N2C">(mol);
    validate_contains<"CN1CCOCC1c2nc(O)c(O)c(n2)C(=O)NCc3ccc(F)cc3">(mol);
    //validate_contains<"CN1CCOC[C@@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2">(mol); // FIXME: stereo
    //validate_contains<"CN1CCOC[C@H]1C2=NC(C(=O)NCc3ccc(F)c(F)c3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    validate_contains<"CN1CCS(=O)CC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
    //validate_contains<"CN1CCS(=O)C[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    //validate_contains<"CN1CC[C@@H](O)C[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    //validate_contains<"CN1C[C@@H](F)C[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    //validate_contains<"CN1C[C@H](CNC(C)=O)OC1=O">(mol); // FIXME: stereo
    //validate_contains<"CN1C[C@H](C[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C)NS(C)(=O)=O">(mol); // FIXME: stereo
    //validate_contains<"CN1C[C@H](F)C[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    validate_contains<"CN1N=C(C(O)=O)C(=O)C2=NC=CC=C12">(mol);
    validate_contains<"CN1N=C(C)N(C)C1=O">(mol);
    validate_contains<"CN1N=C(C2=C1N=CN=C2N)C3=CC=CC=C3">(mol);
    validate_contains<"CN1N=C(C2=CC=C(N)C=C2)C3=C1N=CN=C3N">(mol);
    validate_contains<"CN1N=C(N)C2=CNC3=NC=NC1=C23">(mol);
    validate_contains<"CN1N=CC2=C1C(=O)NCC2">(mol);
    validate_contains<"CN1N=CC2=CC(=CC=C12)C3=CC=CC=C3">(mol);
    validate_contains<"CN1N=NC(C)=C1C">(mol);
    validate_contains<"CN1N=NC2=C1C=C(C=C2)C(N3C=NC=N3)C4=CC=C(Cl)C=C4">(mol);
    validate_contains<"CN1SC(=O)N(CC2=CC=CC=C2)C1=O">(mol);
    //validate_contains<"CN1[C@H]2CC[C@@H]1C[C@H](O)C2">(mol); // FIXME: stereo
    validate_contains<"CN=CCC1=CC=CC=C1">(mol);
    validate_contains<"CN=S">(mol);
    //validate_contains<"CNC(=NS(=O)(=O)C1=CC=C(Cl)C=C1)N2C[C@@H](C(=N2)C3=CC=C(Cl)C=C3)C4=CC=CC=C4">(mol); // FIXME: stereo
    validate_contains<"CNC(=O)C(CC(=O)C(OCC1=CC=CC=C1)C(O)C(O)C(OCC2=CC=CC=C2)C(O)NC3C(C)CC4=C3C=CC=C4)C(C)C">(mol);
    validate_contains<"CNC(=O)C(F)(F)F">(mol);
    validate_contains<"CNC(=O)C1(C)CCC(=O)N1">(mol);
    validate_contains<"CNC(=O)C1=C(C)NC(CC2C(=O)NC3=CC=CC=C23)=C1C">(mol);
    validate_contains<"CNC(=O)C1=CC2=C(C=CC=C2)C=N1">(mol);
    validate_contains<"CNC(=O)C1=CC=C(C)C=C1">(mol);
    validate_contains<"CNC(=O)C1=CC=C(C=C1)C(N)=N">(mol);
    validate_contains<"CNC(=O)C1=CC=C(Cl)C=C1">(mol);
    validate_contains<"CNC(=O)C1=CC=C(Cl)S1">(mol);
}
