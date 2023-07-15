#include "Validate.hpp"

void BindingDB_substructure_part_41(OpenBabel::OBMol &mol)
{
    // SMARTS 2001 - 2050
    validate_contains<"NC(=O)C1=CC=NN1">(mol);
    validate_contains<"NC(=O)C1=CC=[N+](COC[N+]2=C(C=CC=C2)C=NO)C=C1">(mol);
    validate_contains<"NC(=O)C1=CCOC2=CC=CC=C12">(mol);
    //validate_contains<"NC(=O)C1=CNC(C=C2/C(=O)NC3=CC=CC=C23)=C1">(mol); // FIXME: stereo
    validate_contains<"NC(=O)C1=CNC2=C1C=CC=C2">(mol);
    validate_contains<"NC(=O)C1=NC=NC2=C1C(N)=CN2">(mol);
    validate_contains<"NC(=O)C1CC1C(N)=O">(mol);
    validate_contains<"NC(=O)C1CCCC(=O)N1">(mol);
    validate_contains<"NC(=O)C1CNC2CC3=CNC4=CC=CC(C2C1)=C34">(mol);
    validate_contains<"NC(=O)CCC1=CC=CC=C1">(mol);
    validate_contains<"NC(=O)CN1N=NC2=CC=CC=C12">(mol);
    validate_contains<"NC(=O)CNC=O">(mol);
    validate_contains<"NC(=O)N(O)CCC1=CC2=CC=CC=C2S1">(mol);
    validate_contains<"NC(=O)N=NC(N)=O">(mol);
    //validate_contains<"NC(=O)c1ccc[n+](c1)[C@@H]2O[C@H](COP([O-])(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](OP(O)(O)=O)[C@@H]3O)n4cnc5c(N)ncnc45)[C@@H](O)[C@H]2O">(mol); // FIXME: stereo
    validate_contains<"NC(=[NH2+])c1cc2c(I)cccc2s1">(mol);
    validate_contains<"NC(=[NH2+])c1ccc2[nH]c(cc2c1)-c3cccc(-c4ccccc4)c3[O-]">(mol);
    validate_contains<"NC(=[NH2+])c1ccc2[nH]c(nc2c1)-c3ccccc3O">(mol);
    validate_contains<"NC(CC1=CC2=C(N1)C=CC=C2)C(N)=O">(mol);
    validate_contains<"NC(CC1=CC=C(O)C=C1)C(O)=O">(mol);
    validate_contains<"NC(CC1=CC=CC=C1)C(=O)C2=NC3=C(S2)C=CC=C3">(mol);
    validate_contains<"NC(CC1=CN=CN1)C(O)=O">(mol);
    validate_contains<"NC(CCO)C(O)=O">(mol);
    validate_contains<"NC(CO)C=CSC1=CC=CC=C1">(mol);
    validate_contains<"NC(N)(N)S">(mol);
    //validate_contains<"NC(N)=N/C(N)=N">(mol); // FIXME: stereo
    //validate_contains<"NC(N)=N/C1=NCC(CSCCC(N)=NS(N)([O-])[O-])=CS1">(mol); // FIXME: stereo
    validate_contains<"NC1(CCCC1)C2=NC3=NC=NC=C3N=C2">(mol);
    validate_contains<"NC1(CCNCC1)C(O)=O">(mol);
    validate_contains<"NC1=C(C#N)C(C2=CC=CC=C2)=C(C#N)C(=O)N1C3=CC=CC=C3">(mol);
    validate_contains<"NC1=C(C=C(C=C1)N(=O)=O)C(O)=O">(mol);
    validate_contains<"NC1=C(C=C(C=C1)[N+]([O-])=O)C(O)=O">(mol);
    validate_contains<"NC1=C(C=C2C=NC=CC2=N1)C3=C(Cl)C=CC=C3">(mol);
    validate_contains<"NC1=C(C=C2C=NC=CC2=N1)C3=CC=CC=C3">(mol);
    validate_contains<"NC1=C(C=CC=C1)C(O)=O">(mol);
    validate_contains<"NC1=C(C=CC=C1)N(O)O">(mol);
    validate_contains<"NC1=C(C=N)C(C2=CC=CC=C2)=C(C=N)C(=O)N1C3=CC=CC=C3">(mol);
    validate_contains<"NC1=C(N)C=CC=C1">(mol);
    validate_contains<"NC1=C(N=C2C=CC(=CN12)N3CCOCC3)C4=CN=CC=C4">(mol);
    validate_contains<"NC1=C(NC(=O)C2=CC=C(CN3CC(=C)C4=CC=CC=C4C3=O)C=C2)C=CC=C1">(mol);
    validate_contains<"NC1=C2N=C(O)C=CC2=C(OC2=CC=CC=C2)C(O)=C1">(mol);
    validate_contains<"NC1=CC(=NN1)C2=CC=CC=C2">(mol);
    validate_contains<"NC1=CC(=NN1C2=CC=CC=C2)C3CCCC3">(mol);
    validate_contains<"NC1=CC(Br)=CN=C1">(mol);
    validate_contains<"NC1=CC(Br)=CN=C1N">(mol);
    validate_contains<"NC1=CC(Cl)=CC(N)=C1">(mol);
    validate_contains<"NC1=CC(N)=NC=N1">(mol);
    validate_contains<"NC1=CC(NC2=NC=CC(=N2)C2=CC=NC=C2)=CC=C1">(mol);
    validate_contains<"NC1=CC2=C(C=C1)C(=O)NC(N)=N2">(mol);
    validate_contains<"NC1=CC2=C(C=C1)N=CC=C2">(mol);
}
