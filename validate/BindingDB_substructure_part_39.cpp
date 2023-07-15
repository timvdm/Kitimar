#include "Validate.hpp"

void BindingDB_substructure_part_39(OpenBabel::OBMol &mol)
{
    // SMARTS 1901 - 1950
    //validate_contains<"I.N1C=CC=C1">(mol); // FIXME: components
    //validate_contains<"I.NCC1C=NC2=C1C=CC=C2">(mol); // FIXME: components
    validate_contains<"IC1=C(I)C=CC=C1">(mol);
    validate_contains<"IC1=CC(I)=CC=C1">(mol);
    validate_contains<"IC1=CC=CC(CNC2=C3C=CC=CC3=NC4=C2CCCC4)=C1">(mol);
    validate_contains<"IC1=CC=CN1">(mol);
    validate_contains<"N#CC1=C(NC2=CC=CC=C2)C3=CC=CC=C3N=C1">(mol);
    validate_contains<"N(*1****1)C2=NC=NC3=C2C=CC=C3">(mol);
    validate_contains<"N(C1=CC=CC=C1)C1=NC=NC2=C1SC1=C2C=CC=C1">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=CC=NC(NC3=CC=CC=C3)=N2">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC(=CC=N2)C3=CN=CC=C3">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC(=CC=N2)C3=CN=CS3">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC(=CN3C=CN=C23)C4=CC=CC=C4">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC(NC3=CC=CC=C3)=NC=C2">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC=CC(=N2)C3=CN=CS3">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC=CC(OC3=CC=CC=C3)=N2">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC=NC(=N2)C3=CC=CN=C3">(mol);
    validate_contains<"N(C1=NC(=CS1)C1=NC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"N(C1=NC(=CS1)C2=CC=CC=N2)C3=CC=CC=N3">(mol);
    validate_contains<"N(C1=NC=C(S1)C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"N(C1=NC=C(S1)C2=CC=CC=C2)C3=NC=CC=C3">(mol);
    validate_contains<"N(C1=NC=CS1)C2=NC=CC=N2">(mol);
    validate_contains<"N(c1c(cc(cc1)Oc2ccccc2)NC=O)C(=NC(=O)OC)NC(=O)OC">(mol);
    //validate_contains<"N.F.NN.C1=CC2C=CC=CC2C=C1.CCCC(=C)NC3=CCCC3">(mol); // FIXME: components
    //validate_contains<"N.N.[OH-].C1CCCCC1">(mol); // FIXME: components
    //validate_contains<"N.O=N(=O)C1=CC(NC2=NC=CC=N2)=CC=C1">(mol); // FIXME: components
    validate_contains<"N1C(=CC=C1C2=CC=CS2)C3=CC=CS3">(mol);
    validate_contains<"N1C2=C(C=CC=C2)C3=C1C=NC=C3">(mol);
    validate_contains<"N1C2=C(C=CC=C2)C3=CC=NC=C13">(mol);
    validate_contains<"N1C2=CC=CC=C2C3=C1C=C4C=CC=CC4=C3">(mol);
    validate_contains<"N1C2=CC=CC=C2C3=C1C=NC=C3">(mol);
    validate_contains<"N1C2=CC=CC=C2C=CC3=C1C=CC=C3">(mol);
    validate_contains<"N1C2=CC=CC=C2N=CC3=C1C=CC=C3">(mol);
    validate_contains<"N1C2=CC=CN=C2C(=C1C3=CC=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"N1C2=CN=CC=C2C3=C1C=CC=C3">(mol);
    validate_contains<"N1C=C(C2=CC=NN2)C3=C1C=CC(=C3)C4=CC=CC=C4">(mol);
    validate_contains<"N1C=C(C=S1)C2=CC=CC=C2">(mol);
    validate_contains<"N1C=CC(=C1)C2=CC=C3C=CC=CC3=C2">(mol);
    validate_contains<"N1C=CC2=C1C=CC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"N1C=CC2=C1N=CC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"N1C=CC2=C1NC3=C2C=CN3">(mol);
    validate_contains<"N1C=CC2=CC3=CC=CC=C3C=C12">(mol);
    validate_contains<"N1C=CC2=CN=CC=C12">(mol);
    //validate_contains<"N1C=CC=C1.C2=CC=CC=C2.C3=CC=CC=C3">(mol); // FIXME: components
    //validate_contains<"N1C=CC=C1.C2CCCCC2.C3=CC4=C(C=C3)C=CC=C4">(mol); // FIXME: components
    validate_contains<"N1C=CC=C1C2=CC=CC=C2">(mol);
    validate_contains<"N1C=CN=N1">(mol);
    validate_contains<"N1C=NC2=C1C=CC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"N1N=C(C=C1C1=CC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"N1N=CC2=C1C=C3C=CC=CC3=C2">(mol);
}
