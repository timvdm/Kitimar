#include "Validate.hpp"

void BindingDB_substructure_part_33(OpenBabel::OBMol &mol)
{
    // SMARTS 1601 - 1650
    validate_contains<"CNC(=O)C1=CC=C(OCC2CCCCC2)C=C1">(mol);
    validate_contains<"CNC(=O)C1=CN=CS1">(mol);
    validate_contains<"CNC(=O)C1CCC(CC1)C(C)N">(mol);
    validate_contains<"CNC(=O)C1CCCN1C(=O)NC">(mol);
    validate_contains<"CNC(=O)CCC(=O)NC">(mol);
    validate_contains<"CNC(=O)CCCCC1CCSS1">(mol);
    validate_contains<"CNC(=O)CN1C=CC=C1">(mol);
    validate_contains<"CNC(=O)COC(=O)NCC1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CNC(=O)NC">(mol);
    validate_contains<"CNC(=O)NN=C">(mol);
    //validate_contains<"CNC(=O)NN=C1/CCCC2=CC=CN=C12">(mol); // FIXME: stereo
    validate_contains<"CNC(=O)NS(=O)(=O)C1=CC=CC=C1">(mol);
    //validate_contains<"CNC(=O)[C@H](CC(O)=O)NC(C)=O">(mol); // FIXME: stereo
    validate_contains<"CNC(C)C(CC(O)C(N)CC(CC1=CC=C(OC)C(OCCCOC)=C1)C(C)C)C(C)C">(mol);
    validate_contains<"CNC(C)CC1=CC2=C(OCO2)C=C1">(mol);
    validate_contains<"CNC(C)COC1=C(OC)C=C2C(NC3=CC=C(Br)C=C3F)=NC=NC2=C1">(mol);
    validate_contains<"CNC(N)=O">(mol);
    validate_contains<"CNC(N)N">(mol);
    validate_contains<"CNC1=C(C=NC2=CC(F)=C(F)C=C12)C3=NN(C(Cl)=C3C)C4=CC=CC=C4">(mol);
    validate_contains<"CNC1=C(N2CCNCC2)C(Cl)=CC=C1">(mol);
    validate_contains<"CNC1=CC(OC2=CC=C(NC(=O)NC3=CC=CC(=C3)C(F)(F)F)C=C2)=NC=N1">(mol);
    validate_contains<"CNC1=CC(OC2=CC=C(NC(=O)NC3=CC=CC=C3)C=C2)=NC=N1">(mol);
    validate_contains<"CNC1=CC2=C(NC(N)=N2)C=C1C(N)=O">(mol);
    validate_contains<"CNC1=CC=C(C=C1)C(C)=O">(mol);
    validate_contains<"CNC1=CC=C2N(C)C=C(C(N)=O)C(=O)C2=C1">(mol);
    validate_contains<"CNC1=CC=CC(=C1)C2CCCCP2">(mol);
    validate_contains<"CNC1=CC=CC=C1">(mol);
    validate_contains<"CNC1=CC=CC=N1">(mol);
    validate_contains<"CNC1=NC(C)=CS1">(mol);
    validate_contains<"CNC1=NC2=C(C=CC=C2)C(NC)=N1">(mol);
    validate_contains<"CNC1=NC2=C(S1)C(NC3=CC=CC=C3)=NC=N2">(mol);
    validate_contains<"CNC1=NC2=C(SC(=C2)C3=CC=CC(CC(N)C(N)=O)=C3)N4C(C)=CN=C14">(mol);
    validate_contains<"CNC1=NC=CC(=N1)C2=CC=CC=C2">(mol);
    validate_contains<"CNC1=NC=CC(=N1)C2=CC=CN=C2OC3=C(C)C=CC(=C3)C(=O)NC4=CC(=CC=C4N5CCOCC5)C(F)(F)F">(mol);
    validate_contains<"CNC1=NC=CC(NC2=NOC(C)=C2)=N1">(mol);
    validate_contains<"CNC1=NC=CC=N1">(mol);
    validate_contains<"CNC1=NN=C(C)S1">(mol);
    validate_contains<"CNC1=NN=CC=C1">(mol);
    validate_contains<"CNCC(=O)NCC(C)=O">(mol);
    validate_contains<"CNCC(C)=O">(mol);
    validate_contains<"CNCC1=CC=C(N)C=C1">(mol);
    validate_contains<"CNCC1C(=O)N(C)C2=C1C=CC=C2">(mol);
    validate_contains<"CNCCC(OC1=CC=C(C=C1)C(F)(F)F)C2=CC=CC=C2">(mol);
    validate_contains<"CNCCC1=CNC2=C1C(OP(O)(O)=O)=CC=C2">(mol);
    validate_contains<"CNCCCC(NC(=O)C1CCCN1C)C(=O)c2nc3ccccc3s2">(mol);
    validate_contains<"CNS(=C)(=C)C1=C(C=CC=C1)S(=O)(=O)NC">(mol);
    validate_contains<"CNS(=C)(=O)C1=CC=CC2=C1NCC(C)C2">(mol);
    validate_contains<"CNS(=O)(=O)C1=CC=CC2=C1NCC(C)C2">(mol);
    validate_contains<"CNS(=O)(=O)CC1=CC=C2NC=C(CCN(C)C)C2=C1">(mol);
    //validate_contains<"CN[C@H](C(C1=CC=CC=C1)C2=CC=CC=C2)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CCCNC(N)=N)C(=O)C4=NC5=C(S4)C=CC=C5">(mol); // FIXME: stereo
}
