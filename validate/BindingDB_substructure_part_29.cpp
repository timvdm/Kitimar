#include "Validate.hpp"

void BindingDB_substructure_part_29(OpenBabel::OBMol &mol)
{
    // SMARTS 1401 - 1450
    //validate_contains<"CC[C@H]1N(C2CCCC2)C3=NC(NC4=C(OC)C=C(C=C4)C(=O)NC5CCN(C)CC5)=NC=C3N(C)C1=O">(mol); // FIXME: stereo
    validate_contains<"CC[N+](CC)(CC)CCOC1=CC=CC(OCC[N+](CC)(CC)CC)=C1OCC[N+](CC)(CC)CC">(mol);
    validate_contains<"CC[N+](CC)(CC)CCOC1=CC=CC=C1">(mol);
    validate_contains<"CCc1[n+](CC#C)ccn2c(C)ccc12">(mol);
    validate_contains<"CCc1c(C)c2cc3nc(cc4[nH]c(cc5nc(cc1[nH]2)c(C)c5C(C)O)c(C)c4CCC(O)=O)c(CCC(O)=O)c3C">(mol);
    validate_contains<"CCc1nc(N)nc(N)c1-c1ccc2OC(C)(C(=O)N(CCO)c2c1)c1cc(F)cc(F)c1">(mol);
    validate_contains<"CCn1cc(C(=O)C(N)=O)c2ccccc12">(mol);
    validate_contains<"CF">(mol);
    validate_contains<"CN(C(C)=N)C1=CC=CC=C1">(mol);
    validate_contains<"CN(C)C(=O)C(=O)N(C)C(C)(C)C1=NC(C(=O)NCc2ccc(F)cc2)=C(O)C(=O)N1C">(mol);
    validate_contains<"CN(C)C(=O)C(=O)N(C)C1CCN2C1=NC(C(=O)NCc3ccc(F)cc3)=C(O)C2=O">(mol);
    //validate_contains<"CN(C)C(=O)C(=O)N(C)[C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol); // FIXME: stereo
    validate_contains<"CN(C)C(=O)C(=O)NC(C)(C)C1=NC(C(=O)NCc2ccc(F)c(Cl)c2)=C(O)C(=O)N1C">(mol);
    validate_contains<"CN(C)C(=O)C(=O)NC(C)(C)C1=NC(C(=O)NCc2ccc(F)cc2)=C(O)C(=O)N1C">(mol);
    validate_contains<"CN(C)C(=O)C(=O)NC(C)(C)C1=NC(C(=O)NCc2ccc(F)cc2S(C)(=O)=O)=C(O)C(=O)N1C">(mol);
    validate_contains<"CN(C)C(=O)C1=CC(=C([O-])C(=C1)C2=CC=CC=C2)C3=CC4=CC(C(N)=[NH2+])=C(Cl)C=C4N3">(mol);
    validate_contains<"CN(C)C(=O)OC1=CC=CC=C1">(mol);
    validate_contains<"CN(C)C(C)=O">(mol);
    validate_contains<"CN(C)C1=C(C)C=C(N)C=C1">(mol);
    validate_contains<"CN(C)C1=CC=C(C=C1)C(=O)NCC2=CC=C(C=CC(=O)NO)C=C2">(mol);
    validate_contains<"CN(C)C1=CC=C(C=C1)C(=O)NCCCCCC(=O)NO">(mol);
    validate_contains<"CN(C)C1=CC=C(C=C1)C(=O)NCCCCCCC(=O)NO">(mol);
    validate_contains<"CN(C)C1=CC=C(NC(C)=O)C=C1">(mol);
    //validate_contains<"CN(C)C1=CC=C2N=C3C=CC(C=C3SC2=C1)=[N+](/C)C">(mol); // FIXME: stereo
    validate_contains<"CN(C)C1=CC=CC=C1C2=CC=CC=C2">(mol);
    validate_contains<"CN(C)C1=NC=NC2=C1C=CN2C3=CC=CC=C3">(mol);
    validate_contains<"CN(C)C1CC(CO)C(O)C1O">(mol);
    validate_contains<"CN(C)C1CCN(CC2=C(C=C(NC(=O)C3=CC=C(C)C(NC4=NC=CC(=N4)C5=CC=CN=C5)=C3)C=C2)C(F)(F)F)C1">(mol);
    //validate_contains<"CN(C)CC(=O)N1CCN(C)[C@@H](C1)C2=NC(C(=O)NCC3=CC=C(F)C=C3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    //validate_contains<"CN(C)CC(=O)N1CCN(C)[C@@H](C1)C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    validate_contains<"CN(C)CC1=C(C=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CN(C)CC1=CN(N=N1)C2=CC=CC=C2">(mol);
    validate_contains<"CN(C)CC1=CN=NN1C2=C(C)C=NC3=C2C=CC=C3">(mol);
    validate_contains<"CN(C)CCC(CSC1=CC=CC=C1)NC2=CC=C(C=C2[N+]([O-])=O)S(=O)(=O)NC(=O)C3=CC=C(C=C3)N4CCN(CC4)C(=O)C5=C(C=CC=C5)C6=CC=C(Cl)C=C6">(mol);
    validate_contains<"CN(C)CCC(CSC1=CC=CC=C1)NC2=CC=C(C=C2[N+]([O-])=O)S(=O)(=O)NC(=O)C3=CC=C(C=C3)N4CCN(CC5=C(C=CC=C5)C6=CC=C(Cl)C=C6)CC4">(mol);
    validate_contains<"CN(C)CCCC(=O)NC1=CC2=C(C=C1)N=CN=C2NC3=CC=CC=C3">(mol);
    validate_contains<"CN(C)CCCN1C2=CC=CC=C2SC3=C1C=C(C)C=C3">(mol);
    validate_contains<"CN(C)CCCOC1=C2C(C)=CC(C3NC4=C(N3)C=CC=C4)=C2C=CC=C1">(mol);
    validate_contains<"CN(C)CCO">(mol);
    validate_contains<"CN(C)CCOC1=CC=C(C=C1)C(C2=CC(C)=C(O)C(C)=C2)C3=CC(C)=C(O)C(C)=C3">(mol);
    validate_contains<"CN(C)Cc1ccncc1">(mol);
    validate_contains<"CN(C)S(=O)(=O)C1=CC=CC2=C1C=CNC2=O">(mol);
    //validate_contains<"CN(C)S(=O)(=O)CC(=O)N(C)[C@@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol); // FIXME: stereo
    //validate_contains<"CN(C)S(=O)(=O)CC(=O)N(C)[C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol); // FIXME: stereo
    validate_contains<"CN(C)S(=O)(=O)N(C)C1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN(C)S(=O)(=O)N(C)C1CCN2C1=NC(C(=O)NCc3ccc(F)cc3)=C(O)C2=O">(mol);
    //validate_contains<"CN(C)S(=O)(=O)N(C)[C@@H]1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol); // FIXME: stereo
    //validate_contains<"CN(C)S(=O)(=O)N(C)[C@@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol); // FIXME: stereo
    validate_contains<"CN(C)S(=O)(=O)N1CCN(C)C(C1)C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
    validate_contains<"CN(C)S(=O)(=O)NC(C)(C)C1=NC(C(=O)NCc2ccc(F)cc2)=C(O)C(=O)N1C">(mol);
}
