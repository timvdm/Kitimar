#include "Validate.hpp"

void BindingDB_substructure_part_36(OpenBabel::OBMol &mol)
{
    // SMARTS 1751 - 1800
    validate_contains<"COC1=CC=NC=C1C(F)(F)F">(mol);
    validate_contains<"COC1Cc2ncnc(Nc3cccc(O)c3)c2CC1OC">(mol);
    //validate_contains<"COC1O[C@H](CO)[C@@H](O)[C@H](OCC2=CC=CC=C2)[C@H]1OCC3=CC=CC=C3">(mol); // FIXME: stereo
    validate_contains<"COCCNC1=C(C=NC2=CC=CC=C12)C(=O)NNCC3=CC=CC=C3">(mol);
    validate_contains<"COCCNS(=O)(=O)c1ccc(Nc2nccc(n2)-c3c(C)nc4cccnn34)cc1">(mol);
    validate_contains<"COCN1C(C)=NC2=C1C(=O)NC(=O)N2CC(C)C">(mol);
    //validate_contains<"COCOC.CNC(C)CC1=CC=CC=C1">(mol); // FIXME: components
    validate_contains<"COCOC1=C(OCOC)C=C2C(NC3=CC(=CC=C3)C#C)=NC=NC2=C1">(mol);
    validate_contains<"CON=C">(mol);
    validate_contains<"CON=CCC1=CN=CC=C1">(mol);
    validate_contains<"CONC1=CC(O)=CC(NOC)=N1">(mol);
    validate_contains<"COP(=O)(F)OC">(mol);
    validate_contains<"COP(=O)(OC)OC1=CC=CC=C1">(mol);
    validate_contains<"COP(O)(O)=O">(mol);
    //validate_contains<"CO[C@@H]1[C@@H](C[C@H]2O[C@]1(C)[N@@]3C4=CC=CC=C4C5=C6CNC(=O)C6=C7C8=C(C=CC=C8)[N@]2C7=C35)N(C)C(=O)C9=CC=CC=C9">(mol); // FIXME: stereo
    //validate_contains<"COc1cc(ccc1N)C(=O)N[C@@H](C(C)C)C(=O)NC(C)C">(mol); // FIXME: stereo
    //validate_contains<"COc1cc2N(C=C(C(O)=O)C(=O)c2cc1Cc3cccc(Cl)c3F)[C@H](CO)C(C)C">(mol); // FIXME: stereo
    validate_contains<"COc1cc2c(Nc3cnc(NC(=O)c4cccc(Cl)c4)nc3)ncnc2cc1OCCCN5CCOCC5">(mol);
    validate_contains<"COc1cc2ncc3cnc4ccccc4c3c2cc1OC">(mol);
    //validate_contains<"COc1ccc(C=C2/C(=O)N(N=C2C)C(=O)c3ccccc3O)cc1">(mol); // FIXME: stereo
    validate_contains<"COc1ccc(C=C2N(Cc3ccccc3)C(=O)NC2=O)cc1">(mol);
    //validate_contains<"COc1ccc(O)c(C(=O)c2ccc(cc2)C(=O)O[C@@H]3CCC[NH2+]C[C@H]3NC(=O)c4ccncc4)c1F">(mol); // FIXME: stereo
    validate_contains<"COc1ccc(cc1)-c2oc3ncnc(NCCO)c3c2-c4ccc(OC)cc4">(mol);
    validate_contains<"COc1ccc(cc1OCCc2ccc(Cl)cc2Cl)C(=O)NCC3CCN(CC3)C(N)=N">(mol);
    validate_contains<"COc1ccc2[nH]cc(CCNC(=O)CCCCCCNc3c4CCCCc4nc5ccccc35)c2c1">(mol);
    validate_contains<"CP(=O)C1=CC=CC=C1">(mol);
    validate_contains<"CP(C)(S)=O">(mol);
    validate_contains<"CP(S)(=O)C1=CC=C(N)C=C1">(mol);
    validate_contains<"CP(S)(=O)C1=CC=CC=C1">(mol);
    validate_contains<"CS(=O)(=O)C1=CC=C2CN(CC3=CC=CC=C3)CC2=C1">(mol);
    validate_contains<"CS(=O)(=O)CC(=O)NC1CCCN2C(=O)C(O)=C(N=C12)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"CS(=O)(=O)CC[NH2+]Cc1ccc(o1)-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2">(mol);
    validate_contains<"CS(=O)(=O)N=CC1=C(O)C=CC2=C1C=CC=C2">(mol);
    validate_contains<"CS(=O)(=O)NC1=CC=C(O)C2=C1C=CC=C2">(mol);
    //validate_contains<"CS(=O)CC[C@](CO)(C(=O)O[C@H]1CN2CCC1CC2)c3ccccc3">(mol); // FIXME: stereo
    validate_contains<"CS(O)(=O)=O">(mol);
    validate_contains<"CSC">(mol);
    validate_contains<"CSC1=CC=C(C=C1)C2=CC=C(C=C2)S(=O)(=O)N(CC(=O)NO)OC(C)C">(mol);
    validate_contains<"CSC1=NC(C2=CC=C(Cl)C=C2)=C(N=N1)C3=CC=C(Cl)C=C3">(mol);
    validate_contains<"CSC1=NC=CC(=N1)C(O)=O">(mol);
    validate_contains<"CSCC1=CN=CN=C1">(mol);
    validate_contains<"CSCCC(=O)NCCCC(O)=O">(mol);
    validate_contains<"CSCCCC(NC(=O)NC1=CC=CC=C1)C(O)=O">(mol);
    //validate_contains<"CSCC[C@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@@H](NC(=O)[C@@H](N)CS)C(C)C)C(O)=O">(mol); // FIXME: stereo
    //validate_contains<"CSC[C@H](NC(=O)COC1=CC=CC2=C1C=CN=C2)C=O">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H](Cl)CCCl">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H](NC(=O)c1cc(cc(c1)C(=O)N[C@@H](Cc2ccccc2)[C@H](O)CNC3CC3)N(C)S(C)(=O)=O)c4ccccc4">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H](NC(=O)c1cc(cc(c1)C(=O)OCC(N)(CO)Cc2ccccc2)N(C)S(C)(=O)=O)c3ccc(F)cc3">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H](O)CC(C)[C@H](CCC#C)C(C)(C)C">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H](OCC1=CC=CC=C1)C(=O)N[C@@H]2[C@H](O)CC3=CC=CC=C23">(mol); // FIXME: stereo
}
