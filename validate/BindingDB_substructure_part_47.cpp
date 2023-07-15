#include "Validate.hpp"

void BindingDB_substructure_part_47(OpenBabel::OBMol &mol)
{
    // SMARTS 2301 - 2350
    //validate_contains<"O=C1NC2=C(C=CC=C2)C1=C/C3=CC=CN3">(mol); // FIXME: stereo
    //validate_contains<"O=C1NC2=C(C=CC=C2)C1=C/C3=NC4=CC=CC=C4N3">(mol); // FIXME: stereo
    //validate_contains<"O=C1NC2=C(C=CC=C2)C1=C3/CC4=CC=CC=C4N3">(mol); // FIXME: stereo
    validate_contains<"O=C1NC2=C(C=CC=C2)C1=CC3=CC=CN3">(mol);
    validate_contains<"O=C1NC2=C(C=CC=C2)C1=NNC=S">(mol);
    validate_contains<"O=C1NC2=C(C=CC=C2)C1CC3=CC=CN3">(mol);
    validate_contains<"O=C1NC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"O=C1NC2=C(C=CC=C2)C=C1C3=CSC=N3">(mol);
    validate_contains<"O=C1NC2=C(C=NC(NC3=CC=CC=C3)=N2)C=C1C4=CC=CC=C4">(mol);
    validate_contains<"O=C1NC2=C(N1)C3=C(NC4=C3C=CC=C4)C(=O)C5=C2C6=C(N5)C=CC=C6">(mol);
    validate_contains<"O=C1NC2=C(N1)N=CC=C2">(mol);
    validate_contains<"O=C1NC2=C(NC=N2)C(=O)N1">(mol);
    validate_contains<"O=C1NC2=CC(NC3=CC=CC=N3)=CC=C2C=C1">(mol);
    //validate_contains<"O=C1NC2=CC=C(C=C2C1=O)S(=O)(=O)N3CCC[C@H]3COC4=CC=CN=C4">(mol); // FIXME: stereo
    //validate_contains<"O=C1NC2=CC=CC=C2C1=C3/CC4=C(N3)C=CC=C4">(mol); // FIXME: stereo
    validate_contains<"O=C1NC2=CC=CC=C2C1=CC3=CC=CN3">(mol);
    validate_contains<"O=C1NC2=CC=CN=C2N1C3CCCC3">(mol);
    validate_contains<"O=C1NC2=NC(NC3CCNCC3)=NC=C2C=C1C4=CC=CC=C4">(mol);
    validate_contains<"O=C1NC2C(NC3=CC=CC=C23)C4=C1NC5=C4C=CC=C5">(mol);
    validate_contains<"O=C1NC=C2C=CC=C(C2=C1)S(=O)(=O)N3CCCNCC3">(mol);
    validate_contains<"O=C1NC=CC(CC2=CC=CC=C2)=C1">(mol);
    validate_contains<"O=C1NC=CC2=C1C=CC=C2S(=O)(=O)N3CCCNCC3">(mol);
    validate_contains<"O=C1NC=NC2=C1CCCC2">(mol);
    validate_contains<"O=C1NC=NC2=C1CCNC2">(mol);
    validate_contains<"O=C1NC=NN1">(mol);
    validate_contains<"O=C1NCC(NS(=O)(=O)C2=CC=CC=C2)C1CNCC3CC3">(mol);
    validate_contains<"O=C1NCC2=C1C3=C(NC4=C3C=CC=C4)C5=C2C6=C(N5)C=CC=C6">(mol);
    validate_contains<"O=C1NCCC2=C1NN=C2">(mol);
    validate_contains<"O=C1NCCCC2C=CN=C12">(mol);
    validate_contains<"O=C1NCCN1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1NCCO1">(mol);
    validate_contains<"O=C1NCOCN1">(mol);
    validate_contains<"O=C1OC2=C(C=C1)C=C3C=CC=CC3=C2">(mol);
    validate_contains<"O=C1OC2=CC=CC=C2C=C1C3=CC=CC=C3">(mol);
    validate_contains<"O=C1OC=CC2=C1C(=O)OC=C2">(mol);
    validate_contains<"O=C1OCC(=C1C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"O=C1OCC=CO1">(mol);
    validate_contains<"O=C1OCCCCCCC2=CC=CC=C12">(mol);
    validate_contains<"O=C1OCCN1C2=CC=C(C=C2)N3CCOCC3">(mol);
    validate_contains<"O=C1[N]C2=NC=NC=C2N1">(mol);
    validate_contains<"O=CC1=CC2=C(C=CC=C2)N1CC3=CC=CC=C3">(mol);
    validate_contains<"O=CC1=CC=C(C=C1)C2=CC=NC=C2">(mol);
    validate_contains<"O=CC1=CC=CC2=CC=CN=C12">(mol);
    validate_contains<"O=CC1=COC(=N1)C2=CC=CC=C2">(mol);
    validate_contains<"O=CC1CCC2CC3=CNC4=CC=CC(=C34)C2=C1">(mol);
    validate_contains<"O=CC1CCC2CCC3C4CCC5CCCCC5C4C(=O)C=C3C2C1">(mol);
    validate_contains<"O=CCC1=CC=CC=C1">(mol);
    validate_contains<"O=CN1CCCC1C">(mol);
    validate_contains<"O=CN1CCCC1C#N">(mol);
    validate_contains<"O=CNC1=CC=CC=N1">(mol);
}
