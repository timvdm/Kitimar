#include "Validate.hpp"

void BindingDB_substructure_part_59(OpenBabel::OBMol &mol)
{
    // SMARTS 2901 - 2950
    validate_contains<"CCOc1cc(Nc2c(cnc3cc(OCC4CCN(C)CC4)c(OC)cc23)C#N)c(Cl)cc1Cl">(mol);
    //validate_contains<"CC[C@H](N)C(=O)N[C@@H](C(C)C)C(=O)N1CCC[C@H]1C(=O)NC(c2ccccc2)c3ccccc3">(mol); // FIXME: stereo
    validate_contains<"CCc1ccc(cc1)C(=O)Nc2ncc(Nc3ncnc4cc(OCCCN5CCOCC5)c(OC)cc34)cn2">(mol);
    validate_contains<"CN(=O)=O">(mol);
    validate_contains<"CN(C)C(=O)C(=O)N(C)C1CCCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN(C)C(C)(C)c1nc(O)c(O)c(n1)C(=O)NCc2ccc(F)cc2">(mol);
    validate_contains<"CN(C)CC1=NC2=C(SC3=CC=C(C)C=C23)C(=O)N1">(mol);
    validate_contains<"CN(C)CC1=NC2=C(SC3=CC=CC=C23)C(=O)N1">(mol);
    validate_contains<"CN(C)CCOC(C)=O">(mol);
    validate_contains<"CN(C)CCOC(C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    //validate_contains<"CN(C)CC[C@H](CSc1ccccc1)Nc2ccc(cc2[N+]([O-])=O)S(=O)(=O)NC(=O)c3ccc(cc3)N4CCN(CC4)Cc5ccccc5-c6ccccc6">(mol); // FIXME: stereo
    //validate_contains<"CN(C)S(=O)(=O)N(C)[C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol); // FIXME: stereo
    validate_contains<"CN(C)c1ccc(F)c(CCNC(=O)Nc2ccc(Br)cn2)c1F">(mol);
    //validate_contains<"CN(C)c1ccc(cc1)C=C2/C(=O)Nc3ccccc23">(mol); // FIXME: stereo
    //validate_contains<"CN(C1CCCCC1)C(=O)CC[C@@H]([C@H]2CCCCC2)N3Cc4cc(Oc5ccccc5)ccc4N=C3N">(mol); // FIXME: stereo
    validate_contains<"CN(C1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)S(=O)(=O)N4CCN(C)CC4">(mol);
    //validate_contains<"CN(CCCCCCCCN(C)CCC(=O)N1CCC[C@@H]2C3CC4=C(C=CC(=O)N4)[C@]12CC(C)=C3)CCC(=O)N5CCC[C@@H]6C7CC8=C(C=CC(=O)N8)[C@]56CC(C)=C7">(mol); // FIXME: stereo
    validate_contains<"CN(CCCNC(=O)CCc1c[nH]c2ccccc12)CCCNc3c4CCCCc4nc5cc(Cl)ccc35">(mol);
    validate_contains<"CN(Cc1cnc2nc(N)nc(N)c2n1)c3ccc4cc(Cl)ccc4c3">(mol);
    //validate_contains<"CN([C@H]1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3)C(=O)CS(C)(=O)=O">(mol); // FIXME: stereo
    validate_contains<"CN1C(=O)C(O)=C(N=C1C2COCCN2)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CN1C(=O)CCCC1=O">(mol);
    validate_contains<"CN1C(=O)N(C)C2=CC=CC=C2C1=O">(mol);
    validate_contains<"CN1C(=O)N(C)c2cc(ccc12)-c3[nH]c(nc3-c4ccc(F)cc4)-c5ccc(cc5)S(C)=O">(mol);
    validate_contains<"CN1C(=O)N(C)c2nc[nH]c2C1=O">(mol);
    validate_contains<"CN1C(=O)N(C)c2ncn(C)c2C1=O">(mol);
    validate_contains<"CN1C(=O)N(C2=C1C=NC3=CC=C(C=C23)C4=CC5=CC=CC=C5N=C4)C(C)(C)C#N">(mol);
    validate_contains<"CN1C(=O)c2cccnc2N(C3CC3)c4nc(Cl)ccc14">(mol);
    validate_contains<"CN1C(=S)SC(=C)C1=O">(mol);
    validate_contains<"CN1C(C=O)=CC2=C1C=CC=C2">(mol);
    //validate_contains<"CN1C(CC[C@@H](O)C[C@@H](O)CC(O)=O)=C(C(=C1C(=O)NC2=CC(C)=CC=C2)C3=CC=C(F)C=C3)C4=CC=C(F)C=C4">(mol); // FIXME: stereo
    validate_contains<"CN1C=NC2=C(N)N=CN=C12">(mol);
    validate_contains<"CN1C=NC=N1">(mol);
    validate_contains<"CN1CCC(CCOC2=CC=CC=C2)C(O)C1">(mol);
    validate_contains<"CN1CCC=C(*)C1">(mol);
    validate_contains<"CN1CCCC1=O">(mol);
    validate_contains<"CN1CCCC1C2=CC=CC=C2">(mol);
    //validate_contains<"CN1CCCN([C@@H]2CCCCN3C(=O)C(O)=C(N=C23)C(=O)NCc4ccc(F)cc4)S1(=O)=O">(mol); // FIXME: stereo
    validate_contains<"CN1CCN(CC1)c2ccc(Nc3ncc4C=C(C(=O)N(C)c4n3)c5c(Cl)cccc5Cl)cc2">(mol);
    validate_contains<"CN1CCN(CC1)c2ccc(cc2)C(=O)Nc3n[nH]c4cn(cc34)C(=O)Cc5cccs5">(mol);
    validate_contains<"CN1CCOCC1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol);
    //validate_contains<"CN1CCOC[C@H]1C2=NC(C(=O)NCc3ccc(F)cc3)=C(O)C(=O)N2C">(mol); // FIXME: stereo
    validate_contains<"CN1CCOc2cc3NC(=O)C=C(c3cc12)C(F)(F)F">(mol);
    validate_contains<"CN1N=CC2=CC=CC=C12">(mol);
    //validate_contains<"CN=C(/N=C/NC1=CC=CC=C1)C2=CC=CC=C2">(mol); // FIXME: stereo
    validate_contains<"CNC(=O)C(=O)CCCCCCC(=O)Nc1cccc(c1)-c2ccccc2">(mol);
    //validate_contains<"CNC(=O)C1=C(C)NC(C=C2/C(=O)NC3=CC=CC=C23)=C1C">(mol); // FIXME: stereo
    validate_contains<"CNC(=O)C1=CC=CS1">(mol);
    validate_contains<"CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN">(mol);
    //validate_contains<"CNC(=O)OCCCCCN1[C@H](Cc2ccccc2)[C@H](O)[C@@H](O)[C@@H](Cc3ccccc3)N(CCCCCOC(=O)NC)C1=O">(mol); // FIXME: stereo
}
