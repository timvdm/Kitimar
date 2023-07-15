#include "Validate.hpp"

void BindingDB_substructure_part_69(OpenBabel::OBMol &mol)
{
    // SMARTS 3401 - 3450
    validate_contains<"N1N=CC2=CC=CC=C12">(mol);
    validate_contains<"NC1=C2N=CNC2=NC=N1">(mol);
    validate_contains<"NC1=CC=NC(N)=N1">(mol);
    validate_contains<"O=C1NC2=NC=NC=C2C=C1">(mol);
    validate_contains<"O=CNC1=NC=CC=C1">(mol);
    validate_contains<"OCCCC1=CNC2=CC=CC=C12">(mol);
    validate_contains<"N1C=CC=N1">(mol);
    validate_contains<"ONC">(mol);
    validate_contains<"C1=CN(C=N1)C2=CC=CC=C2">(mol);
    //validate_contains<"C1CCCC1.C2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"CBr">(mol);
    validate_contains<"CN(C(C)=O)C1=CC=CC=C1">(mol);
    validate_contains<"CN1CCN(CC1)Cc2ccc(cc2)C(=O)Nc3ccc(C)c(Nc4nccc(n4)-c5cccnc5)c3">(mol);
    validate_contains<"C[S]">(mol);
    validate_contains<"NC1=CC=CN=C1N">(mol);
    validate_contains<"NC1CCCCC1">(mol);
    //validate_contains<"O=C1NC2=CC=CC=C2\C1=C\C3=CC=CN3">(mol); // FIXME: stereo
    validate_contains<"O=C1OCCN1C2=CC=CC=C2">(mol);
    validate_contains<"OC(=O)C1=CC=CC=C1">(mol);
    validate_contains<"OC1CCCNC1">(mol);
    validate_contains<"CCNCCC1=CC=CC=C1">(mol);
    validate_contains<"O=C1NCCCN1">(mol);
    validate_contains<"C1=CC2=C(C=C1)C=NC=C2">(mol);
    validate_contains<"C1=CC2=C(C=C1)N=CC=C2">(mol);
    validate_contains<"N=C=S">(mol);
    validate_contains<"*NC">(mol);
    validate_contains<"C1CC2=CC=CC=C2O1">(mol);
    validate_contains<"C1COC2=CC=CC=C2C1">(mol);
    validate_contains<"CC(=O)NC1=CC=CC=C1">(mol);
    validate_contains<"CC12CCC3C(CCC4=C3C=CC(O)=C4)C1CCC2O">(mol);
    validate_contains<"CCCC1=CC=CC=C1">(mol);
    validate_contains<"CCl">(mol);
    validate_contains<"N1N=CC2=C1C=CC=C2">(mol);
    validate_contains<"SC1=NC=CC=C1">(mol);
    validate_contains<"CNCC">(mol);
    validate_contains<"O=C1NCCCCN1">(mol);
    validate_contains<"C1C2=CC=CC=C2C3=C1C=NN3">(mol);
    validate_contains<"C1CC1">(mol);
    validate_contains<"O1C=CC2=CC=CC=C12">(mol);
    validate_contains<"OCC">(mol);
    validate_contains<"[H]N">(mol);
    validate_contains<"C1CCCCC1">(mol);
    validate_contains<"C1CC2=CC=CC=C2C1">(mol);
    validate_contains<"C1=CN=CN=C1">(mol);
    validate_contains<"CN1C=CC=N1">(mol);
    validate_contains<"O=C1NC2=CC=CN=C2N1">(mol);
    validate_contains<"N1C=NC2=C1C=CC=C2">(mol);
    validate_contains<"O=[N+](O)c1ccccc1">(mol);
    validate_contains<"CCCC">(mol);
    //validate_contains<"C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O">(mol); // FIXME: stereo
}
