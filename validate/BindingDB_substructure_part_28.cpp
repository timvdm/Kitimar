#include "Validate.hpp"

void BindingDB_substructure_part_28(OpenBabel::OBMol &mol)
{
    // SMARTS 1351 - 1400
    validate_contains<"CCN(CC)Cc1c(O)ccc(C=NNC(N)=S)c1O">(mol);
    validate_contains<"CCN(CC1=CN=C2C=CC=CC2=C1)C(=O)CC3=CC=CC=C3">(mol);
    validate_contains<"CCN1C(=NC2=C1C=C(NC(=O)C1=CC=C(OCCN3CCOCC3)C=C1)N=C2)C1=NON=C1N">(mol);
    validate_contains<"CCN1C(=O)CC(C)(C)N=C1N">(mol);
    validate_contains<"CCN1C=C(C(=C1)C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"CCN1C=C(C(=O)C(C)=O)C2=C1C=CC=C2">(mol);
    validate_contains<"CCN1C=C(CC)C2=C1C=CC=C2">(mol);
    validate_contains<"CCN1C=C(NC)C2=C1C=CC=C2">(mol);
    validate_contains<"CCN1CC(CCN2CCOCC2)C(C1=O)(C3=CC=CC=C3)C4=CC=CC=C4">(mol);
    validate_contains<"CCN1CC2=CC=CC=C2CC1C">(mol);
    validate_contains<"CCN1CCCC(O)C1">(mol);
    validate_contains<"CCN1CCCC1CNC(=O)C2=CC(=CC=C2OC)S(N)(=O)=O">(mol);
    validate_contains<"CCN1CCN(CC1)C2=NC(=O)C(=C(N2)C)CC3=CC=CC=C3Cl">(mol);
    validate_contains<"CCN1N=CC=C1C2=CNC3=C2C=C(C=N3)C4=CC(C(=O)N(C)C)=C(N)C=C4">(mol);
    validate_contains<"CCNC(=O)C(CC(O)=O)NC(C)=O">(mol);
    validate_contains<"CCNC(=O)C1=C(OC(C)=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CCNC(=O)C1OC(C(O)C1O)N2C=NC3=C(N)N=CN=C23">(mol);
    validate_contains<"CCNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC">(mol);
    validate_contains<"CCNC1=CC=CC=C1">(mol);
    validate_contains<"CCNCC">(mol);
    validate_contains<"CCO">(mol);
    validate_contains<"CCOC">(mol);
    validate_contains<"CCOC(=O)C(CCC1=CC=CC=C1)NC(C)C(=O)N2CCCC2C(O)=O">(mol);
    validate_contains<"CCOC(=O)C(O)CC1=CC2=C(OCO2)C=C1[N+]([O-])=O">(mol);
    validate_contains<"CCOC(=O)C1=C(COCCN)NC(C)=C(C1C2=CC=CC=C2Cl)C(=O)OC">(mol);
    validate_contains<"CCOC(=O)C1=C(Cl)C2=CC(C)=C(C)C=C2N=C1">(mol);
    validate_contains<"CCOC(=O)C1=C(O)C2=C(C=CC=C2)N=C1">(mol);
    validate_contains<"CCOC(=O)C1=CC(OC(CC)CC)C(NC(C)=O)C(N)C1">(mol);
    validate_contains<"CCOC(=O)C1=CC2=C(S1)C=CC(NC3=NN=CC4=CC(OCCCC5CCCCC5)=C(OC)C=C34)=C2">(mol);
    validate_contains<"CCOC(=O)C1=CNC=C(C1C2=CC=CC=C2)C(=O)OCC">(mol);
    validate_contains<"CCOC(=O)C1=NC(=NO1)C2=CC=C(O)C=C2">(mol);
    validate_contains<"CCOC(=O)C1=NOC(=C1)C2=C(O)C=CC=C2">(mol);
    validate_contains<"CCOC(=O)C1=NOC(=C1)C2=CC=C(C(=C2)N(=O)=O)C(C)(C)C">(mol);
    validate_contains<"CCOC(=O)C1=NOC(=C1)C2=CC=C(O)C=C2">(mol);
    validate_contains<"CCOC(=O)C1=NOC(=C1)C2=CC=C(OC(C)(C)C)C(=C2)N(=O)=O">(mol);
    //validate_contains<"CCOC(=O)C=C[C@H](C[C@@H]1CCNC1=O)NC(=O)[C@H](CC#C)N2C=CC=C(NC(=O)C3=NOC(C)=C3)C2=O">(mol); // FIXME: stereo
    validate_contains<"CCOC(=O)CC(=C)C(O)C1=C(C[N+]([O-])=O)C=CC=C1">(mol);
    validate_contains<"CCOC1=C(OCC)C=C(C=C1)C2=NC(=CS2)C3=NC=CC=C3">(mol);
    //validate_contains<"CCOC1=C(OCCN[C@H](C)CC2=CC(=C(OC)C=C2)S(N)(=O)=O)C=CC=C1">(mol); // FIXME: stereo
    validate_contains<"CCOC1=CC=C(C=C1)C(N)=N">(mol);
    validate_contains<"CCOC1=CNC=C1">(mol);
    validate_contains<"CCOC1=NC2=C(N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NNN=N5)C(=CC=C2)C(O)=O">(mol);
    validate_contains<"CCOCC1C2CCC(CC1C3=CC=C(Cl)C(Cl)=C3)N2C">(mol);
    validate_contains<"CCOCN1C2=NC=CC=C2NC(=O)C3=C1N=CC=C3">(mol);
    validate_contains<"CCOCOCC">(mol);
    validate_contains<"CCSC1=CN=C(NC(C)=O)S1">(mol);
    //validate_contains<"CC[C@@H](C)Cl">(mol); // FIXME: stereo
    //validate_contains<"CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](CC(O)=O)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CCC(O)=O)NC(=O)CNC(=O)[C@H](C)N)[C@@H](C)CC)C(C)C)[C@@H](C)CC)C(=O)N[C@@H](C)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CCC(N)=O)C(N)=O">(mol); // FIXME: stereo
    //validate_contains<"CC[C@H](CO)N1C=C(C(O)=O)C(=O)c2cc(Cc3cccc(Cl)c3F)ccc12">(mol); // FIXME: stereo
    //validate_contains<"CC[C@H](NC)C(=O)N[C@H]1CCC[C@H]2CC[C@H](N2C1=O)C(=O)NC(c3ccccc3)c4ccccc4">(mol); // FIXME: stereo
}
