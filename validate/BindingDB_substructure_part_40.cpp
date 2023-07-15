#include "Validate.hpp"

void BindingDB_substructure_part_40(OpenBabel::OBMol &mol)
{
    // SMARTS 1951 - 2000
    validate_contains<"N1N=CC2=C1C=CC3=C2NN=C3">(mol);
    validate_contains<"N1N=CC2=C1C=NC=N2">(mol);
    validate_contains<"N1N=CC2=C1N=CC=C2">(mol);
    validate_contains<"N1N=NC2=CC=C(C=C12)C3=C(N=NN3)C4=CC=CC=N4">(mol);
    validate_contains<"N1N=NC=C1C2=CC=CC=N2">(mol);
    validate_contains<"N1NC(C2=CNC=C2)C(=C1)C3=CC=CC=C3">(mol);
    //validate_contains<"N1OC(=N/C2=CC=CC=C2)C3=C1C=CC=C3">(mol); // FIXME: stereo
    validate_contains<"N=C1C(=O)NC2=C1C=CC=C2">(mol);
    validate_contains<"N=C1NC(=N)C2=C1C=CC=C2">(mol);
    validate_contains<"N=C1NCC2=C1C=CC=C2">(mol);
    validate_contains<"N=C1NCCN1">(mol);
    validate_contains<"N=S">(mol);
    validate_contains<"NC(*)C(C=C)C1=CC=CC=C1">(mol);
    validate_contains<"NC(=C)C1=C(S)C2=C(N1)C=CC=C2">(mol);
    validate_contains<"NC(=C)C1=C(S)C2=C(N1)C=CC=N2">(mol);
    validate_contains<"NC(=N)C1=CC(CN2C(=CC3=C2C=CC=C3O)C(=O)NCC4=CC=CC=C4)=CC=C1">(mol);
    validate_contains<"NC(=N)C1=CC2=C(NC(=N2)C3=CC=CC(C4=CC=CC=C4)=C3O)C=C1F">(mol);
    validate_contains<"NC(=N)C1=CC2=C(NC=N2)C=C1">(mol);
    validate_contains<"NC(=N)C1=CC=C(C=C1)C2=CC=C(O2)C3=CC=C(C=C3)C(N)=N">(mol);
    validate_contains<"NC(=N)C1=CC=C(C=CC2=CC=C(C=C2)C(N)=N)C=C1">(mol);
    validate_contains<"NC(=N)C1=CC=C(N)C=C1">(mol);
    validate_contains<"NC(=N)C1=CC=C(OCCCCCOC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"NC(=N)C1=CC=C2NC(=CC2=C1)C3=CC=CC(OC4CCCCC4)=C3O">(mol);
    validate_contains<"NC(=N)C1=CC=C2NC(CC2=C1)C3=CC=CC=C3O">(mol);
    validate_contains<"NC(=N)C1=CC=CC(=C1)N2N=C(C=C2C(N)=O)C(F)(F)F">(mol);
    validate_contains<"NC(=N)C1=CC=CC(CN2CCCCC2CCO)=C1">(mol);
    //validate_contains<"NC(=N)N1CC=C(CCNC(=O)[C@@H]2C[C@@H]3CC[C@@H](O)C[C@@H]3N2C(=O)[C@@H](Cc4ccccc4)NC(=O)[C@H](O)Cc5ccccc5)C1">(mol); // FIXME: stereo
    validate_contains<"NC(=N)N1CCCCC1">(mol);
    validate_contains<"NC(=N)NC1=CC=CC=C1">(mol);
    validate_contains<"NC(=N)SCCC1=CC(CCSC(N)=N)=CC=C1">(mol);
    validate_contains<"NC(=N)c1ccc2nc([nH]c2c1)-c3cc4ccccc4cn3">(mol);
    //validate_contains<"NC(=N/C1=NC(=CC2=CC=CC=C12)C3=NC=CC=C3)C4=CC=CC=C4">(mol); // FIXME: stereo
    validate_contains<"NC(=O)C1=C(N)C2=CC=CN=C2S1">(mol);
    validate_contains<"NC(=O)C1=C(NC2=C(Cl)C=NC(NC3=CC=C(C=C3)N4CCNCC4)=N2)C=CC=C1">(mol);
    validate_contains<"NC(=O)C1=C(NCC2=C3C=CC=CC3=NC=C2)C=CS1">(mol);
    validate_contains<"NC(=O)C1=C(S)C2=C(N1)C=CC=N2">(mol);
    validate_contains<"NC(=O)C1=C(SC2=CC=CC=C2)C3=C(N1)C=CC=N3">(mol);
    validate_contains<"NC(=O)C1=CC(=NN1)C(F)(F)F">(mol);
    validate_contains<"NC(=O)C1=CC(O)=NC=C1">(mol);
    validate_contains<"NC(=O)C1=CC2=C(N1)C=CC=C2">(mol);
    validate_contains<"NC(=O)C1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"NC(=O)C1=CC=C(C=C1)C2=NC=C(N2)C3=NC=CC=C3">(mol);
    validate_contains<"NC(=O)C1=CC=C(Cl)C=C1">(mol);
    validate_contains<"NC(=O)C1=CC=C(Cl)S1">(mol);
    validate_contains<"NC(=O)C1=CC=C(N)C=C1">(mol);
    validate_contains<"NC(=O)C1=CC=C(OC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"NC(=O)C1=CC=CC2=C1N=CN=C2">(mol);
    validate_contains<"NC(=O)C1=CC=CC2=C1NC3=CC=CC=C3C2=O">(mol);
    validate_contains<"NC(=O)C1=CC=CN1">(mol);
    validate_contains<"NC(=O)C1=CC=CN=C1">(mol);
}
