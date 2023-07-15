#include "Validate.hpp"

void BindingDB_substructure_part_30(OpenBabel::OBMol &mol)
{
    // SMARTS 1451 - 1500
    //validate_contains<"CN(C)[C@@H]1CCN2C1=NC(C(=O)NCc3ccc(F)cc3)=C(O)C2=O">(mol); // FIXME: stereo
    validate_contains<"CN(C)c1ccc(cc1)-c2[nH]nc3-c4cccc(NC(=O)CN5CCOCC5)c4C(=O)c23">(mol);
    //validate_contains<"CN(C)c1cccc2c(cccc12)S(=O)(=O)N[C@@H](Cc3cccc(c3)[CH+](N)N)C(=O)N4CCN(CC4)S(C)(=O)=O">(mol); // FIXME: stereo
    validate_contains<"CN(C1=CC(Cl)=CC=C1)S(=O)(=O)C2=CC3=C(NC(=O)C3CC4=C(C)C(=C(C)N4)C(=O)N5CCNCC5)C=C2">(mol);
    validate_contains<"CN(C1=CC=C(C=C1)C(C(F)(F)F)C(F)(F)F)S(=O)(=O)C2=CC=CC=C2">(mol);
    //validate_contains<"CN(C1CCCC1)C=C/C=C2CNCN2">(mol); // FIXME: stereo
    validate_contains<"CN(C1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)C(=O)CS(C)(=O)=O">(mol);
    validate_contains<"CN(C1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)C(C)=O">(mol);
    validate_contains<"CN(C1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(C)(=O)=O">(mol);
    validate_contains<"CN(C1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)C(=O)C(F)F">(mol);
    validate_contains<"CN(C1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(=O)(=O)N4CCOCC4">(mol);
    validate_contains<"CN(C1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(C)(=O)=O">(mol);
    validate_contains<"CN(C1CCN2C1=NC(C(=O)NCc3ccc(F)cc3)=C(O)C2=O)C(C)=O">(mol);
    validate_contains<"CN(CC(C)=O)CC(C)=O">(mol);
    //validate_contains<"CN(CC1=CC=CC=C1)[C@H]2CCCN3C(=O)C(O)=C(N=C23)C(=O)NCC4=CC=C(F)C=C4">(mol); // FIXME: stereo
    validate_contains<"CN(CCCC1=CC=C(CC2CC(=O)CC2=O)C=C1)C3=CC=CC=C3">(mol);
    validate_contains<"CN(CCOC1=CC=C(CC2SC(=O)NC2=O)C=C1)C3=NC=CC=C3">(mol);
    //validate_contains<"CN(Cc1ccccc1)[C@H]2CCCN3C(=O)C(O)=C(N=C23)C(=O)NCc4ccc(F)cc4">(mol); // FIXME: stereo
    //validate_contains<"CN([C@@H]1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(C)(=O)=O">(mol); // FIXME: stereo
    //validate_contains<"CN([C@@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(=O)(=O)N4CCCC4">(mol); // FIXME: stereo
    //validate_contains<"CN([C@@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(C)(=O)=O">(mol); // FIXME: stereo
    //validate_contains<"CN([C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCC3=CC=C(F)C=C3)C(C)=O">(mol); // FIXME: stereo
    //validate_contains<"CN([C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)C(=O)Cn4cccn4">(mol); // FIXME: stereo
    //validate_contains<"CN([C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)C(C)=O">(mol); // FIXME: stereo
    //validate_contains<"CN([C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(=O)(=O)N4CCN(C)CC4">(mol); // FIXME: stereo
    //validate_contains<"CN([C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(=O)(=O)c4c(C)nn(C)c4C">(mol); // FIXME: stereo
    validate_contains<"CN1C(=CC2=C1C=CC=C2)C(N)=O">(mol);
    validate_contains<"CN1C(=O)C(=C(C)C2=CC=CC=C12)C(=O)NCC(O)=O">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)N2CCOCC2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)C(=O)N2CCOCC2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2c[nH]cn2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2ccccn2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2cccnn2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2cscn2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2ncccn2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2nnc(C)[nH]2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(=O)c2nnc(C)o2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NC(C)=O)C(=O)NCc2ccc(F)cc2">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C(C)(C)NS(C)(=O)=O)C(=O)NCc2ccc(F)cc2">(mol);
    validate_contains<"CN1C(=O)C(O)=C(N=C1C2COCCN2C(=O)OC(C)(C)C)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)C2(NN=C(S2)C3=CC=CC=C3)C4=CC=CC=C14">(mol);
    validate_contains<"CN1C(=O)C2=CC=CC=C2N=C1C=CC1=CC=CC=C1">(mol);
    validate_contains<"CN1C(=O)C2N(CCC3=C2NC4=C3C=CC=C4)C1=O">(mol);
    //validate_contains<"CN1C(=O)C=C(N2CCC[C@@H](N)C2)N(Cc3ccccc3C#N)C1=O">(mol); // FIXME: stereo
    validate_contains<"CN1C(=O)C=C(O)C2=C1C=CC=C2">(mol);
    validate_contains<"CN1C(=O)C=CC2=CN=C(NC3=CC=CC=C3)N=C12">(mol);
    validate_contains<"CN1C(=O)CCN=C1N">(mol);
    validate_contains<"CN1C(=O)COCC1=O">(mol);
    validate_contains<"CN1C(=O)N(C)C2=CC=CC=C12">(mol);
    //validate_contains<"CN1C(=O)N(C)[C@@]2(Oc3c(O)cc4c(CC[N+]4(C)C)c3C5=C2C(=O)c6c(CCN)c[nH]c6C5=O)C1=O">(mol); // FIXME: stereo
}
