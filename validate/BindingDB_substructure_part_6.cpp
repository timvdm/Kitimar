#include "Validate.hpp"

void BindingDB_substructure_part_6(OpenBabel::OBMol &mol)
{
    // SMARTS 251 - 300
    validate_contains<"CNC1=C2C=C(OC)C=CC2=NC2=C1C=CC=C2">(mol);
    validate_contains<"CNC1=C2NC=NC2=NC(NC2=CC=CC=C2)=N1">(mol);
    validate_contains<"CNC1=CC(=NC(N)=N1)C1=CC=CC=C1">(mol);
    validate_contains<"CNC1=CC=CN2C=NN=C12">(mol);
    validate_contains<"CNC1=NC(C)=CC2=NC=CN12">(mol);
    validate_contains<"CNCC12CC3CC(C1)CC(C3)(C2)OC">(mol);
    validate_contains<"CNS(=O)(=O)C1=CC=CC2=CN=CC=C12">(mol);
    validate_contains<"COC(=O)C1=CC(=CC(C)=C1OC)C(=CCCN)C1=CC(C(=O)OC)=C(OC)C(C)=C1">(mol);
    validate_contains<"COC(=O)C1=CC=C(C=C1)N=C=O">(mol);
    validate_contains<"COC(=O)CC(C)N">(mol);
    validate_contains<"COC(=O)NC1=NC2=CC=CC=C2N1">(mol);
    validate_contains<"COC1=C(O)C=C(O)C2=C1OC(CC2=C)C1=CC=CC=C1">(mol);
    validate_contains<"COC1=C(O)C=C(O)C2=C1OC(CC2=O)C1=CC=CC=C1">(mol);
    validate_contains<"COC1=C(O)C=CC(=C1)C1OCC2C1COC2C1=CC(OC)=C(O)C=C1">(mol);
    validate_contains<"COC1=C(O)C=CC(C=CCO)=C1">(mol);
    validate_contains<"COC1=C(OC)C=C(C=C1)C1=NNC(=O)C2CCCCC12">(mol);
    validate_contains<"COC1=C(OC)C=C2C(=C1)N=CN=C2N1CCN(CC1)C(=O)OCC1=CC=CC=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(NC3=C2C(=CNC3=O)C2=CC=C(O)C=C2)=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(OC3=CC=C(N)C=C3)=CC=NC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(OC3=CC=C(NC(N)=O)C=C3)=CC=NC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(OC3=CC=C(NC(N)=S)C=C3)=CC=NC2=C1">(mol);
    validate_contains<"COC1=C2C=CC(=O)NC2=CC=C1">(mol);
    validate_contains<"COC1=C2C=CC=NC2=CC=C1">(mol);
    //validate_contains<"COC1=CC(=CC(OC)=C1OC)[C@H]1[C@@H]2[C@H](COC2=O)C(O)C2=CC3=C(OCO3)C=C12">(mol); // FIXME: stereo
    validate_contains<"COC1=CC(C=CC(=O)OCC2=CC=CC=C2)=CC=C1O">(mol);
    validate_contains<"COC1=CC(NC2=C(C=NC3=C2C=C(OC)C(OCCCN2CCN(C)CC2)=C3)C#N)=C(Cl)C=C1Cl">(mol);
    validate_contains<"COC1=CC2=C(C=C1OC)C(N)=NC(=N2)N1CCNCC1">(mol);
    validate_contains<"COC1=CC=C(C=C1)C1=NNC(=O)C1=NNC1=CC=CC=C1">(mol);
    validate_contains<"COC1=CC=C(C=C1)S(=O)(=O)N1CCN(CC1)C(C)C1=NC2=C(C=CC=C2)C(OC2CCCCC2)=N1">(mol);
    validate_contains<"COC1=CC=C(C=C1O)C1=CC2=NC=CN2C(NC2=NC=CC=C2C(N)=O)=N1">(mol);
    validate_contains<"COC1=CC=C(C=C1OC)C1=CC2=NC=CN2C(NC2=NC=CC=C2C(N)=O)=N1">(mol);
    validate_contains<"COC1=CC=C(CCN2CCC(CC2)NC2=NC3=C(C=CC=C3)N2CC2=CC=C(F)C=C2)C=C1">(mol);
    validate_contains<"COC1=CC=CC2=C1C(=O)C1=C(O)C3=C(CC(O)CC3OC3CC(N)C(O)C(C)O3)C(O)=C1C2=O">(mol);
    validate_contains<"COC1=NC2=C(N)C=C(O)C(OC3=CC=CC=C3)=C2C(C)=C1">(mol);
    validate_contains<"COC1=NC2=C(NC(C)CCCN)C=C(O)C(OC3=CC(=CC=C3)C(F)(F)F)=C2C(C)=C1">(mol);
    validate_contains<"COC1=NC2=C(NC(C)CCCN)C=C(OC)C(OC3=CC(=CC=C3)C(F)(F)F)=C2C(C)=C1">(mol);
    validate_contains<"COCCOC1=C(OC)C=C2C(=C1)N=CC1=C2C=CN=C1N">(mol);
    validate_contains<"COCCOC1=CC2=C(C=C1OCCOC)C(NC1=CC(=CC=C1)C#C)=NC=N2">(mol);
    validate_contains<"CON1C(=O)C(=CC2=CN=C(N)N=C12)C1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"CON1C(=O)C(=CC2=CN=C(NC3=CC=C(C=C3)N3CCN(C)CC3)N=C12)C1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"COc1cc(N)ccc1-c1ccc2c(n[nH]c2c1)-c1nc2c(cccc2[nH]1)N1CCN(C)CC1">(mol);
    validate_contains<"CS(=O)(=O)CCNCc1ccc(o1)-c1ccc2ncnc(Nc3ccc(OCc4cccc(F)c4)c(Cl)c3)c2c1">(mol);
    validate_contains<"CS(=O)(=O)N1CCC2=C(C1)C=CS2">(mol);
    validate_contains<"CS(=O)(=O)NC1=CC=CC2=C1N=CC=C2">(mol);
    validate_contains<"CSC1=NC2=C(C=CC=C2)C2=C1N=CN2">(mol);
    validate_contains<"CSC1=NC=C(C(=O)N2CCC(F)C2)C(=N1)N1CCC(N)CC1">(mol);
    validate_contains<"CSCCC(NC(=O)NC1=CC=CC=C1)C(O)=O">(mol);
    //validate_contains<"C[As+](C)(C)CCO[C@@H](C1=CC=CC=C1[N+](O)=O)C(F)(F)F">(mol); // FIXME: stereo
    //validate_contains<"C[C@@H](CNC(C)(C)CC1=CC=C(OC2=NC=CC=C2)C=C1)COC3=CC=CC4=C3C5=C(N4)C=CC=C5">(mol); // FIXME: stereo
    //validate_contains<"C[C@H]1CC=C[C@@H](O)[C@@H](O)CC=C/CC2=C(C(O)=CC(O)=C2)C(=O)O1">(mol); // FIXME: stereo
}
