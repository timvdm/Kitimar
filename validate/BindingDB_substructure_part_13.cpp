#include "Validate.hpp"

void BindingDB_substructure_part_13(OpenBabel::OBMol &mol)
{
    // SMARTS 601 - 650
    validate_contains<"C1=CC2=C(C=C1)C3=C4C5=CC=CC=C5C6=C(C7=CC=CC=C7C=C6)C4=C8C=C9C=CC=CC=CC=C(C=CC=CC=CC=C)C=CC=CC9=CC8=C3C=C2">(mol);
    validate_contains<"C1=CC2=C(C=C1)C=C(C=C2)C3=CC4=C(C=CC=C4)C=C3">(mol);
    //validate_contains<"C1=CC2=C(C=C1)C=CC=C2.C3=CC4=C(C=C3)C=CC=C4.C5=CC6=C(C=C5)C=CC=C6">(mol); // FIXME: components
    validate_contains<"C1=CC2=C(C=CN=C2)C=N1">(mol);
    validate_contains<"C1=CC2=C(C=CN=C2)N=C1">(mol);
    validate_contains<"C1=CC2=C(C=N1)C=NC=C2">(mol);
    validate_contains<"C1=CC2=C3C(C=CC4=CC=CC(C=C2)=C34)=C1">(mol);
    validate_contains<"C1=CC2=CC3=CC4=CC=CC=C4C=C3C=C2C=C1">(mol);
    validate_contains<"C1=CC2=CC3=CC=CC=C3C=C2C=C1">(mol);
    validate_contains<"C1=CC2=CC=C3C=CC4=C5C3=C2C6=C7C(C=CC(C=C4)=C57)=CC=C16">(mol);
    validate_contains<"C1=CC2=CC=C3C=CC=C4C=CC(=C1)C2=C34">(mol);
    validate_contains<"C1=CC2=NC=CN2N=C1">(mol);
    validate_contains<"C1=CC=C(C=C1)C#CC2=CC=C3N=CC=CC3=C2">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=C3C=CC=CC3=NC=N2">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=CC(=NC(=C2)C3=CC=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=CC=C3C=CC=CC3=C2">(mol);
    validate_contains<"C1=CC=C(C=C1)C2=CN=C3C=CN=CC3=C2">(mol);
    validate_contains<"C1=CC=C(C=C1)N2C=CC=C2">(mol);
    validate_contains<"C1=CC=C(C=C1)N2C=CC=C2C3=CC=CC=C3">(mol);
    validate_contains<"C1=CC=C(C=C1)N=NC2=CC=CC=C2">(mol);
    validate_contains<"C1=CC=C2C(=C1)C=CC3=CC=CC=C23">(mol);
    validate_contains<"C1=CC=C2C(=C1)C=CC3=NC=CN=C23">(mol);
    validate_contains<"C1=CC=C2C(C=CC3=C2C4=CC=CC=C4C=C3)=C1">(mol);
    validate_contains<"C1=CC=C2C(C=CC3=CC=CC=C23)=C1">(mol);
    //validate_contains<"C1=CC=CC=C1.C2=CC3=CC=CC=C3C=C2">(mol); // FIXME: components
    //validate_contains<"C1=CC=CC=C1.C2=CC=CC=C2">(mol); // FIXME: components
    //validate_contains<"C1=CC=CC=C1.C2=CC=CC=C2.C3=CC=CC=C3">(mol); // FIXME: components
    //validate_contains<"C1=CC=CC=C1.C2=CC=CC=C2.C3=CC=CC=C3.C4=CC=CC=C4.C5=CC=CC=C5">(mol); // FIXME: components
    //validate_contains<"C1=CC=CC=C1.CNC(=O)NC2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"C1=CC=NN=C1">(mol);
    validate_contains<"C1=CN(C=C1)C2=CC=CC3=C2C=CC=C3">(mol);
    validate_contains<"C1=CN(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"C1=CN(N=C1)C2=NC(=CC=N2)C3=CC=CN=C3">(mol);
    validate_contains<"C1=CN2C(C=NC3=C2C=CC=C3)=N1">(mol);
    validate_contains<"C1=CN2C=CN=C2C=C1">(mol);
    validate_contains<"C1=CN2C=CN=C2N=C1">(mol);
    validate_contains<"C1=CN2N=CC=C2N=N1">(mol);
    validate_contains<"C1=CN=C2C(=C1)C=CC3=CC=CN=C23">(mol);
    validate_contains<"C1=CN=NN=C1">(mol);
    validate_contains<"C1=NC2=CN=CN=C2N=C1">(mol);
    validate_contains<"C1=NN=C2C=CC=CC12">(mol);
    validate_contains<"C1=NN=CN=N1">(mol);
    validate_contains<"C1=NN=NN=C1">(mol);
    validate_contains<"C1C(=NOC12CCNCC2)C3=CC=CC=C3">(mol);
    validate_contains<"C1C2=C(SC3=C1C=CC=C3)C=CC=C2">(mol);
    validate_contains<"C1C2=CC=CC=C2NC3=C1C=CC=C3">(mol);
    validate_contains<"C1C2=CNN=C2C3=CC=CC=C13">(mol);
    validate_contains<"C1C2C(NN=C2C3=C1C=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"C1C2CC3CC1CC(C2)C3">(mol);
    validate_contains<"C1C2CCCCC2C3C1C4=C5C3CCC6CCC7=CC=CC(=C4)C7=C56">(mol);
}
