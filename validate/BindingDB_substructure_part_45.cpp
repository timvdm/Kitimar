#include "Validate.hpp"

void BindingDB_substructure_part_45(OpenBabel::OBMol &mol)
{
    // SMARTS 2201 - 2250
    //validate_contains<"O=C(C1=C(OC)C=C(O)C(C/C=C(C)C)=C1O)/C=C/C2=CC=C(O)C=C2">(mol); // FIXME: stereo
    validate_contains<"O=C(C1=CN(C2=C1C=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"O=C(C=CC1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"O=C(CC1=CC=CC=C1)N(CCCC2=CC=CC=C2)CC3=CC=CC=C3">(mol);
    validate_contains<"O=C(CC1=CC=CC=C1)NC2=NNC=C2">(mol);
    validate_contains<"O=C(CNC1CCNCC1)N2CCSC2">(mol);
    validate_contains<"O=C(CSC1=NN=CN1C2=C3C=CC=CC3=CC=C2)NC4=CC=CC=C4">(mol);
    validate_contains<"O=C(Cc1ccccc1)NC(=S)Nc2ccccc2">(mol);
    //validate_contains<"O=C(N([H])C([H])([H])[H])[C@]2(N([H])C(C([H])([H])[H])=O)C([H])([H])c1c(c([H])c([H])c([H])c1[H])C2([H])[H]">(mol); // FIXME: stereo
    validate_contains<"O=C(N1C=CC=N1)C2=CC=CC=C2">(mol);
    validate_contains<"O=C(N1CCCCC1)C2=CC3=C(C=C2)C=C(C4=CSC(=N4)C5=CC=NC=C5)C(=O)N3">(mol);
    validate_contains<"O=C(N1CCN(CC1)C2C3=CC=CC=C3C4=C2C=CC=C4)C5=CC6=C(NC=C6)C=C5">(mol);
    //validate_contains<"O=C(N2CC[C@]([H])(CN1CCC[C@]([H])(C1)C(=O)N(CC)CC)CC2)C=5C3=CC=CC=C3C=C4C=CC=CC4=5">(mol); // FIXME: stereo
    validate_contains<"O=C(NC(=O)C1=NN=CC=C1)C2=CC=CN=N2">(mol);
    validate_contains<"O=C(NC(C1=CC=CC=C1)C2=CC=CC=C2)C3CCCN3">(mol);
    validate_contains<"O=C(NC1=**=**=*1)C2=*N=*S2">(mol);
    validate_contains<"O=C(NC1=CC=C(C=C1)C2=CC=CC=C2)C3=CC4=C(O3)C=CC=C4">(mol);
    validate_contains<"O=C(NC1=CC=C(OC2=CC=CC=C2)C=C1)C3=C(NCC4=CC=NC=C4)N=CC=C3">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C(=O)NC2=CC=CC=C2">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=C(NCC3=CC=NC=C3)C=CS2">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CC3=CC=CC=C3C=C2">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CC=CC=C2NCC3=C4C=CC=CC4=NC=C3">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CC=CC=C2SCC3=CC=NC=C3">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CN=C(NC3=NC=NC=C3)S2">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)C2=CN=CS2">(mol);
    validate_contains<"O=C(NC1=CC=CC=C1)NC2=CC=CC=C2">(mol);
    validate_contains<"O=C(NC1=NC(=CS1)C1=NC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"O=C(NC1CC1)NC2=CNN=C2C3=NC4=C(N3)C=CC(CN5CCOCC5)=C4">(mol);
    validate_contains<"O=C(NC1CCCC1)C2=CC3=C(C=CC=C3)C=N2">(mol);
    validate_contains<"O=C(NCC1=CC=CC=C1)C2=CC=CN2">(mol);
    validate_contains<"O=C(Nc1ccccc1)c1n[nH]c2ccccc12">(mol);
    validate_contains<"O=C(Nc1ccccc1)c2cc[nH]n2">(mol);
    validate_contains<"O=C1C(=O)C2=C(C=CC(=C2)N(=O)=O)C3=C1C=CC=C3">(mol);
    validate_contains<"O=C1C(=O)C2=C(C=CC=C2)C3=CC=CC=C13">(mol);
    validate_contains<"O=C1C2=C(NC3=C1C=CC=C3)C=CC=C2">(mol);
    validate_contains<"O=C1C=C2C=CC=CC2C=N1">(mol);
    validate_contains<"O=C1C=C2C=CC=NC2C=N1">(mol);
    validate_contains<"O=C1C=C2C=CC=NN2C=N1">(mol);
    validate_contains<"O=C1C=CC(=O)C2=CC=CC=C12">(mol);
    validate_contains<"O=C1C=CC(=O)C=C1">(mol);
    validate_contains<"O=C1C=CC1=O">(mol);
    validate_contains<"O=C1C=CC=CN1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1C=CCN1C2=CC=CC=C2">(mol);
    validate_contains<"O=C1C=CNCN1">(mol);
    validate_contains<"O=C1C=COC2=CC=CC=C12">(mol);
    validate_contains<"O=C1CC(=O)NC(=O)N1">(mol);
    validate_contains<"O=C1CC(=O)NC2=CC=CC=C2N1">(mol);
    validate_contains<"O=C1CC(C(=O)N1)C2=CC=CC=C2">(mol);
    validate_contains<"O=C1CC(OC2=C1C=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"O=C1CC2=C(NC3=C2C=CC=C3)C4=CC=CC=C4N1">(mol);
}
