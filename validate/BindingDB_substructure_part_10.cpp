#include "Validate.hpp"

void BindingDB_substructure_part_10(OpenBabel::OBMol &mol)
{
    // SMARTS 451 - 500
    //validate_contains<"O=C1NC2=C(C=CC=C2)C1=N/NC1=CC=CC=C1">(mol); // FIXME: stereo
    validate_contains<"O=C1NC2=C(C=CC=C2)C1=NNC1=CC=CC=C1">(mol);
    validate_contains<"O=C1NC2=C(NC3=C1C=CN=N3)C=CC=C2">(mol);
    validate_contains<"O=C1NC2=C(OC3=C1C=CN=N3)C=CC=C2">(mol);
    validate_contains<"O=C1NC2=CC=C3N=CSC3=C2C1=CNC1=CC=CC=C1">(mol);
    //validate_contains<"O=C1NC2=CC=C3N=CSC3=C2\C1=C\NC1=CC=CC=C1">(mol); // FIXME: stereo
    validate_contains<"O=C1NC2=CC=CC(OC3CCCCC3)=C2C=C1">(mol);
    validate_contains<"O=C1NC=C(C2=C1NC1=CC=CC=C21)C1=CC=CC=C1">(mol);
    validate_contains<"O=C1NC=CC2=C1C=CC=C2">(mol);
    //validate_contains<"O=C1NN=C(C1=N/NC1=CC=CC=C1)C1=CC=CC=C1">(mol); // FIXME: stereo
    validate_contains<"O=C1OC2=C(C=CC=C2)C(=C1)N1C=CC=C1">(mol);
    validate_contains<"O=C1OC2=CC3=C(C=C2C=C1)C1=C(CCCC1)O3">(mol);
    validate_contains<"O=C=NC1=CC=CC=C1">(mol);
    validate_contains<"O=CC1CC(=O)C2=CC=CC=C12">(mol);
    validate_contains<"O=CNC1=CC(NC2=NC=CC(=N2)C2=CC=CN=C2)=CC=C1">(mol);
    validate_contains<"O=CNC1=CC2=C(C=NN2)C=C1">(mol);
    validate_contains<"O=CNC1=NNC2=C1CN(C2)C=O">(mol);
    validate_contains<"O=N(=O)C1=CN2CCCOC2=N1">(mol);
    validate_contains<"O=S(=O)(NC1=CC(NC2=NC=CC(=N2)C2=CC=NC=C2)=CC=C1)C1=C2C=CN=CC2=CC=C1">(mol);
    validate_contains<"O=S(=O)(NC1=CC=NN1)C1=CC=CC2=CN=CC=C12">(mol);
    validate_contains<"O=S(=O)NC1=CC(NC2=NC=CC(=N2)C2=CC=NC=C2)=CC=C1">(mol);
    validate_contains<"OC(=O)C(=O)C(Br)Br">(mol);
    validate_contains<"OC(=O)C(=O)C(Cl)Cl">(mol);
    validate_contains<"OC(=O)C(=O)C(F)(F)F">(mol);
    validate_contains<"OC(=O)C(NC=O)C1=CC=CC=C1">(mol);
    validate_contains<"OC(=O)C1=C(O)C=CC(NCC2=CC(O)=CC=C2O)=C1">(mol);
    validate_contains<"OC(=O)C1CC(=O)C2=C1C=CC=C2">(mol);
    validate_contains<"OC(=O)CCCNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1">(mol);
    validate_contains<"OC(=O)CCCNC(=O)CCSCC1=CC=CC=C1">(mol);
    //validate_contains<"OC(=O)CCCNC(=O)[C@@H](CSCC1=CC=CC=C1)NC=O">(mol); // FIXME: stereo
    validate_contains<"OC(=O)CCP(CCC(O)=O)CCC(O)=O">(mol);
    validate_contains<"OC(=O)CN1C(=O)N(CCCC2=CC=CC=C2)C2=C(C3=C(CCC3)S2)C1=O">(mol);
    validate_contains<"OC(=O)NC1=NNC2=CC(=CC(C3=NC4=CC=CC=C4N3)=C12)C1=CC=CO1">(mol);
    validate_contains<"OC(CN1CCC(O)(CC1)C1=CC(F)=CC=C1)C1=CC2=C(C=C1)NC(=O)CC2">(mol);
    validate_contains<"OC1=C(C=C(Br)C=C1)C1=CC(C2=CC=CC=C2)=C(C#N)C(=O)N1">(mol);
    validate_contains<"OC1=C(Cl)C=C(Cl)C=C1SC1=C(O)C(Cl)=CC(Cl)=C1">(mol);
    validate_contains<"OC1=C(O)C=C(C=C(C#N)C#N)C=C1">(mol);
    validate_contains<"OC1=C2C=CC(=O)OC2=CC=C1">(mol);
    validate_contains<"OC1=CC(O)=C2C(=O)C(O)=C(OC2=C1)C1=CC(O)=C(O)C=C1">(mol);
    validate_contains<"OC1=CC(O)=C2C(=O)C=COC2=C1">(mol);
    validate_contains<"OC1=CC2=C(NCC2)C=C1">(mol);
    //validate_contains<"OC1=CC=C(C=C1O)C(=O)C[N@@+]12CC[C@@H](CC1)[C@@H](C2)NC(=O)C1=CC=CC(Cl)=C1">(mol); // FIXME: stereo
    validate_contains<"OC1=CC=C(CC2=NN=C3N=CC(=NN23)C2=CC=CC=C2)C=C1">(mol);
    validate_contains<"OC1=CC=C2C(=O)C(O)=C([*])OC2=C1">(mol);
    validate_contains<"OC1=CC=C2C(=O)C(O)=C([*])OC2=C1[*]">(mol);
    validate_contains<"OC1=CC=C2C(=O)C=C([*])OC2=C1[*]">(mol);
    validate_contains<"OC1=CC=C2NC(=O)CC2=C1">(mol);
    validate_contains<"OC1=CC=NC=C1C(F)(F)F">(mol);
    validate_contains<"OC1=NC2=CC=C(O)C(OC3=CC=CC=C3)=C2C=C1">(mol);
    validate_contains<"OC1=NC2=CC=CC(OC3=CC=CC=C3)=C2C=C1">(mol);
}
