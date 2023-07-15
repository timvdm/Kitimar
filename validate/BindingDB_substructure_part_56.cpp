#include "Validate.hpp"

void BindingDB_substructure_part_56(OpenBabel::OBMol &mol)
{
    // SMARTS 2751 - 2800
    validate_contains<"c1cc2ccccc2[nH]1">(mol);
    validate_contains<"c1ccc2cc3cc4ccccc4cc3cc2c1">(mol);
    validate_contains<"c1ccccc1NOCCC">(mol);
    validate_contains<"c1ccccn1">(mol);
    validate_contains<"c1cccnc1C(=O)O">(mol);
    validate_contains<"BrC1=CC(NC2=NC=NC3=C2C=CC=C3)=CC=C1">(mol);
    validate_contains<"C(C1CCCCN1)C1=CC=CC=C1">(mol);
    validate_contains<"C(C1CCCN1)C1=CC=CC=C1">(mol);
    validate_contains<"C(C1NCCC2=CC=CC=C12)N1CCCC1">(mol);
    validate_contains<"C1=C(N=C2N=CC=CN12)C1=CC=CC=C1">(mol);
    validate_contains<"C1=CC=C(C=C1)C1=CC2=C(C=CC=C2)C(=N1)C1=CC=CC=C1">(mol);
    validate_contains<"C1=CC=C(C=C1)C1=CC=CC2=C1N=CC=C2">(mol);
    validate_contains<"C1=NN2C=CC=NC2=C1">(mol);
    validate_contains<"CC(=O)NC1=CC(NC2=NC=CC(=N2)C2=CC=CN=C2)=CC=C1">(mol);
    validate_contains<"CC1=CC2=C(C=C1N)N=C(N)N2">(mol);
    validate_contains<"CC1=CC=C(C=C1)N1N=C(C=C1NC(=O)NC1=CC=CC2=CC=CC=C12)C(C)(C)C">(mol);
    validate_contains<"CC1=CC=CC=C1">(mol);
    validate_contains<"CCC1=CNC=C1">(mol);
    validate_contains<"CN1CCCCC1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(=C1)N=CN=C2N1CCNCC1">(mol);
    validate_contains<"COC1=CC2=C(C=C1OCCCN1CCOCC1)C(NC1=CC=C(F)C(Cl)=C1)=NC=N2">(mol);
    validate_contains<"C[N+]1=CC=CC=C1">(mol);
    validate_contains<"NC(CC1=CNC2=C1C=CC=C2)C(O)=O">(mol);
    validate_contains<"NC1=CC(=N)C2=C(C=CC=C2)C1=O">(mol);
    validate_contains<"NC1=CC=NC=C1">(mol);
    validate_contains<"NC1CCNC1">(mol);
    //validate_contains<"N[C@@H](CCC(=O)NC(CS)C(=O)NCC(O)=O)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"O1N=CC2=C1C=CC=C2">(mol);
    validate_contains<"O1N=CC2=CC=CC=C12">(mol);
    validate_contains<"O=C1CCC2=C1C=CC=C2">(mol);
    validate_contains<"O=C1CNC(=O)N1">(mol);
    validate_contains<"OB1OCC2=CC=CC=C12">(mol);
    //validate_contains<"Br.C1=CC=CC=C1">(mol); // FIXME: components
    validate_contains<"BrC1=CC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"C#C">(mol);
    validate_contains<"C(C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"C(OC1=CC=CC=C1)C2CCCN2SC3=CC=C4N(CC5=CC=CC=C5)CCC4=C3">(mol);
    validate_contains<"C1(=NN=C(S1)CCC2=CC=CC=C2)NC(=O)C3=CC=CC=C3">(mol);
    validate_contains<"C1*CC2=C1C=CC3=C2C4=CC=CC=C4*3">(mol);
    validate_contains<"C1=CC2=CC3=C(C=CC=C3)C=C2C=C1">(mol);
    validate_contains<"C1=CC2=CC=C3C=CC=CC3=C2C=C1">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=CC=C3">(mol);
    validate_contains<"C1=CC=C(C=C1)N2N=CC=C2C3=CC=CC=C3">(mol);
    validate_contains<"C1=CC=C2C=CC=CC2=C1">(mol);
    //validate_contains<"C1=CC=CC=C1.CC2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"C1C2=CC=CC=C2C3=C1C=CC=C3">(mol);
    validate_contains<"C1CC1C2=CC=CC3=C2C=CC=C3">(mol);
    validate_contains<"C1CC2=C(C1)C3=C(C=C2)C4=CC=CC=C4C=C3">(mol);
    validate_contains<"C1CC2=C(C1)C=CC=C2">(mol);
    validate_contains<"C1CC2=C(NC3=CC=CC=C13)C=CC=C2">(mol);
}
