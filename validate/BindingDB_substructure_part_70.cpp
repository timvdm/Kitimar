#include "Validate.hpp"

void BindingDB_substructure_part_70(OpenBabel::OBMol &mol)
{
    // SMARTS 3451 - 3488
    //validate_contains<"N1C=CC=C1.C2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"N1C=NC2=CC=CC=C12">(mol);
    validate_contains<"O=C">(mol);
    validate_contains<"OC(=O)C1=C(O)C=CC=C1">(mol);
    validate_contains<"NC1=NC=CC=N1">(mol);
    validate_contains<"O=C1NC2=CC=CC=C2N1">(mol);
    validate_contains<"N1C=CC2=CC=CC=C12">(mol);
    validate_contains<"S1C=NC2=CC=CC=C12">(mol);
    validate_contains<"CCC">(mol);
    validate_contains<"FC">(mol);
    validate_contains<"CN1C=CC2=C1C=CC=C2">(mol);
    validate_contains<"N">(mol);
    validate_contains<"C1CNCCN1">(mol);
    validate_contains<"C1CCC2=CC=CC=C2C1">(mol);
    validate_contains<"C1CCNC1">(mol);
    validate_contains<"N1C=CC2=C1C=CC=C2">(mol);
    validate_contains<"C">(mol);
    validate_contains<"BrC1=CC=CC=C1">(mol);
    validate_contains<"CI">(mol);
    validate_contains<"OC">(mol);
    validate_contains<"S1C2=CC=CC=C2N=C1C1=CC=CC=C1">(mol);
    validate_contains<"CN">(mol);
    validate_contains<"O=S">(mol);
    validate_contains<"[H]C(=C([H])C1=CC(O)=CC(O)=C1)C2=CC=C(O)C=C2">(mol);
    //validate_contains<"N1C=CC=C1.C2=CC3=C(C=C2)C=CC=C3">(mol); // FIXME: components
    validate_contains<"[H]CN">(mol);
    validate_contains<"C1NCC2=C1C=CC3=C2C4=CC=CC=C4N3">(mol);
    validate_contains<"N1C=CC=C1">(mol);
    validate_contains<"OC1=CC=C(C=CC2=CC(O)=CC(O)=C2)C=C1">(mol);
    validate_contains<"C1=CC=NC=C1">(mol);
    validate_contains<"CC12CCC3C(CCC4=CC(O)=CC=C34)C1CCC2O">(mol);
    validate_contains<"CNC">(mol);
    validate_contains<"NC">(mol);
    validate_contains<"C1=CC=CC=C1">(mol);
    validate_contains<"C1=CC2=C(C=C1)C=CC=C2">(mol);
    validate_contains<"O=C1OC2=CC=CC=C2C=C1">(mol);
    validate_contains<"CC">(mol);
    validate_contains<"CC(C)C1(C)SC(Nc2ccccc2C(F)(F)F)=NC1=O">(mol);
}
