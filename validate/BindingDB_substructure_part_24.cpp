#include "Validate.hpp"

void BindingDB_substructure_part_24(OpenBabel::OBMol &mol)
{
    // SMARTS 1151 - 1200
    validate_contains<"CC1=NC(CN)=C(N)O1">(mol);
    validate_contains<"CC1=NC2=C(C)C=CC=C2N1">(mol);
    validate_contains<"CC1=NC2=C(C=C(N3CCOCC3)C(F)=C2)C(N)=C1">(mol);
    validate_contains<"CC1=NC2=C(C=CC=C2)C(N)=C1">(mol);
    validate_contains<"CC1=NC2=C(C=CC=C2)C=N1">(mol);
    validate_contains<"CC1=NC2=C(C=NC=C2)C=C1C3=CC=CC=C3">(mol);
    validate_contains<"CC1=NC2=C(CN1)C=CC=C2">(mol);
    validate_contains<"CC1=NC2=C(N)N=CC=C2O1">(mol);
    validate_contains<"CC1=NC2=C(N1)C=CC(=C2)C3=CC=CC=C3">(mol);
    validate_contains<"CC1=NC2=C(N1)C=CC=C2C">(mol);
    validate_contains<"CC1=NC2=C(OC3=C2C=CC=N3)C=N1">(mol);
    validate_contains<"CC1=NC2=CC(C)=CC=C2N1">(mol);
    validate_contains<"CC1=NC2=CC=CC=C2CS1">(mol);
    validate_contains<"CC1=NC2=CC=CC=C2N1">(mol);
    validate_contains<"CC1=NC=CC2=C1C(C)(C)NC3=C2C=CC(Br)=C3">(mol);
    validate_contains<"CC1=NC=NC2=C1C=CC=C2">(mol);
    //validate_contains<"CC1=NN(C(=O)C1=C/C2=CC=CC=C2)C3=CC=CC=C3">(mol); // FIXME: stereo
    validate_contains<"CC1=NN(C(N)=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CC1=NN(C=C1)C2=CC=CC=C2">(mol);
    //validate_contains<"CC1=NN=C(C[C@H](C[C@H](O)CN2CCN(CC3=CC=C(O3)C4=CC=C(Cl)C=C4)C[C@H]2C(=O)NCC(F)(F)F)C(=O)N[C@@H]5[C@H](O)COC6=C5C=CC=C6)O1">(mol); // FIXME: stereo
    validate_contains<"CC1=NN=C(N)C2=CC=CC=C12">(mol);
    validate_contains<"CC1=NN=C(O1)C1=C(C=CS1)S(=O)(=O)C1=CC=CC=C1">(mol);
    validate_contains<"CC1=NN=C(O1)C1=CC=CS1">(mol);
    validate_contains<"CC1=NN=C(O1)C2=CC=C(C)C(C)=C2">(mol);
    validate_contains<"CC1=NNC(=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CC1=NNC(C)=C1C(O)=O">(mol);
    validate_contains<"CC1=NOC2(C1)CCNCC2">(mol);
    validate_contains<"CC1=NOC2=C1C(=O)NN=C2CC3=CC=CC=C3">(mol);
    validate_contains<"CC1C(NC2=CC(Cl)=NC(SCC(O)=O)=N2)C=CC=C1C">(mol);
    //validate_contains<"CC1C(O)CCC2(C)C1CCC3(C)C2C(O)CC4C(C(CC34C)OC(C)=O)=C(CCC=C(/C)C)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"CC1C(S)N(C)C2=C1C=CC=C2">(mol);
    validate_contains<"CC1C2=CC=CC=C2SC3=C1C=CC=C3">(mol);
    validate_contains<"CC1C2=CC=CC=C2SCC3=C1C=CC=C3">(mol);
    validate_contains<"CC1C2CC3=CC=CC=C3C2=NN1C">(mol);
    validate_contains<"CC1C2CCCC2=NN1C">(mol);
    validate_contains<"CC1C2CN(CC12)C3=CC(=C(C=C3)C#N)C(F)(F)F">(mol);
    validate_contains<"CC1C=C(NC(=O)C2=CC(NC(=O)C3=CC(=CC3C)N(C)C(=O)c4ccc(cc4)C(=O)N(C)C5=CC(C)C(=C5)C(=O)NC6=CC(C)C(=C6)C(=O)NC7=CC(C)C(=C7)C(=O)NCCCN(C)C)=CC2C)C=C1C(=O)NCCCN(C)C">(mol);
    validate_contains<"CC1CC(=O)N(C)C1=O">(mol);
    validate_contains<"CC1CC(=O)N(C1=O)C2=CC=C(O)C=C2">(mol);
    validate_contains<"CC1CC(C)C2=C(C1)C(C)C=CC2C">(mol);
    validate_contains<"CC1CC(N(C)C1)C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
    //validate_contains<"CC1CC(O)C=C(C)C2CC(CC12)=C(/C)C">(mol); // FIXME: stereo
    validate_contains<"CC1CC23CC(C)CC24CC(C)CC45CC(C)CC56CC(C)CC36C1">(mol);
    validate_contains<"CC1CC2=C(CN1C)C=C(C)C=C2">(mol);
    validate_contains<"CC1CC2=C(CN1C)C=C3OC(C)C(=O)NC3=C2">(mol);
    validate_contains<"CC1CC2=C(CN1C)C=C3OCC(=O)N(C)C3=C2">(mol);
    validate_contains<"CC1CC2=CC3=C(OC(C)CO3)C=C2CN1">(mol);
    validate_contains<"CC1CC2C(CCC3=CC(O)=CC=C23)C4CCC(O)C14">(mol);
    validate_contains<"CC1CCC(C1)CC2CCCC2">(mol);
    validate_contains<"CC1CCC2(C1)CC(C)CC3(CCCC3)C2">(mol);
}
