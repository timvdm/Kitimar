#include "Validate.hpp"

void BindingDB_substructure_part_16(OpenBabel::OBMol &mol)
{
    // SMARTS 751 - 800
    validate_contains<"C1NC2=C(C=NC(NC3=CC=CC=C3)=N2)C=C1C4=CC=CC=C4">(mol);
    validate_contains<"C1NC2=CC=CC=C2C1C3=C(CNC3)C4CC5=C(N4)C=CC=C5">(mol);
    validate_contains<"C1NC2=CC=CC=C2C1C3=C(CNC3)C4CNC5=C4C=CC=C5">(mol);
    validate_contains<"C1NC2=CC=CC=C2C3=CC=CC=C13">(mol);
    validate_contains<"C1NC2=CC=CC=C2CN3C=CC=C13">(mol);
    validate_contains<"C1NC2=CC=CC=C2N1">(mol);
    validate_contains<"C1NC2=CC=CC=C2N=C1">(mol);
    validate_contains<"C1NC2=CN=NN2C3=C1C=CN3">(mol);
    validate_contains<"C1NCC(=C1)C2=CNC3=C2C=CC=C3">(mol);
    validate_contains<"C1NCC2=C(N1)SC=C2">(mol);
    validate_contains<"C1NCC2=C1C3=C(NC4=C3CCCC4)C5=C2C6=C(N5)C=CC=C6">(mol);
    validate_contains<"C1NCC2=C1C3=C(NC4=CC=CC=C34)C5=NC=CC=C25">(mol);
    validate_contains<"C1NCC2=C1C=CC3=C2C4=C(CCCC4)N3">(mol);
    validate_contains<"C1NCC2=C1CCC3NC4=CC=CC=C4C23">(mol);
    validate_contains<"C1NCC2=C1NN=C2">(mol);
    validate_contains<"C1NNN=C1">(mol);
    validate_contains<"C1N[H]OC[H]1">(mol);
    validate_contains<"C1OC2=C(C=C1)C=CC=C2">(mol);
    validate_contains<"C1OC2=C(O1)C=C(C=C2)C3=CNC=N3">(mol);
    validate_contains<"C1OC2=C(OC1C3=NC=CS3)C=CC=C2">(mol);
    validate_contains<"C1OC2=CC=CC=C2OC1C3=CC=CC=C3">(mol);
    validate_contains<"C1SC2=NN=CN2N=C1">(mol);
    validate_contains<"C2N1COCN1CO2">(mol);
    validate_contains<"C=C">(mol);
    validate_contains<"C=C(C1=CC=CC=C1)C2=CC3=C(NC=C3)C=C2">(mol);
    validate_contains<"C=C(NS(=O)(=O)C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    //validate_contains<"C=C.CCC">(mol); // FIXME: components
    validate_contains<"C=C1C(=O)NC(=O)C2=C1C=CC=C2">(mol);
    validate_contains<"C=C1CNC(C2=CC=CC=C2)C3=CC=CC=C3N1">(mol);
    validate_contains<"C=C1CNC2=CC=CC=C2N1">(mol);
    validate_contains<"C=C1N(SC2=CC=CC=C12)C3=CC=CC=C3">(mol);
    validate_contains<"C=C1NC(=O)CS1">(mol);
    validate_contains<"C=C1SC(=S)NC1=O">(mol);
    validate_contains<"C=CC(=N)NC1=CC=C2N=CN=C(NC3=CC=CC=C3)C2=C1">(mol);
    validate_contains<"C=CC(=O)NC1=CC=C2C=CC=C(NC3=CC=CC=C3)C2=C1">(mol);
    validate_contains<"C=CC(=O)NC1=CC=C2N=CN=C(NC3=CC=CC=C3)C2=C1">(mol);
    validate_contains<"C=CC1=CC2=C(NC=C2)C=C1">(mol);
    validate_contains<"C=CC1=CC=CC=C1C2=CC=CN2">(mol);
    //validate_contains<"C=CC=C/C(C=C)=C/C=C.CC=C/C1=CC=C2C=CC(=CC)C(=C/C)C2=C1C">(mol); // FIXME: components
    //validate_contains<"C=CC=C/C=C.C=CC=C1/C(=C)C=CC2=CC=CC=C12">(mol); // FIXME: components
    validate_contains<"CB(O)O">(mol);
    validate_contains<"CC#CC1(O)CCC2C3CCC4=CC(=O)CCC4=C3C(CC12C)C5=CC=C(C=C5)N(C)C">(mol);
    validate_contains<"CC#CCC(C)C(O)CCC1C(O)CC2CC(CCCCC(O)=O)CC12">(mol);
    validate_contains<"CC(=C)C1=C(C)C=C(C(C)=C1)C2=CC3=C(C=CC=C3)C=C2">(mol);
    validate_contains<"CC(=C)C1=C(S)C2=C(N1)C=CC=C2">(mol);
    validate_contains<"CC(=C)C1=C(S)C2=C(N1)C=CC=N2">(mol);
    validate_contains<"CC(=C)C1=CC2=C(N1)C=CC=C2">(mol);
    validate_contains<"CC(=C)C1=CC2=C(N1)C=CC=N2">(mol);
    validate_contains<"CC(=C)C1CCC(C)=CC1">(mol);
    //validate_contains<"CC(=C/CO)C#CC(O)C(C)=CCC1=C(C)CCCC1(C)C">(mol); // FIXME: stereo
}
