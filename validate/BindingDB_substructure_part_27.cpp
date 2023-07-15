#include "Validate.hpp"

void BindingDB_substructure_part_27(OpenBabel::OBMol &mol)
{
    // SMARTS 1301 - 1350
    validate_contains<"CCCCC(CC)N1[Se]C2=CC=CC=C2C1=O">(mol);
    validate_contains<"CCCCC1=C(CC)C=C(CCC)C=C1C(=O)NC">(mol);
    validate_contains<"CCCCC1=C(CC)C=CC=C1">(mol);
    validate_contains<"CCCCC1=CC(CC)=NC(CO)=C1">(mol);
    validate_contains<"CCCCC1=CC=CC=C1">(mol);
    validate_contains<"CCCCC1=CNC=C1">(mol);
    validate_contains<"CCCCC1=CNC=N1">(mol);
    validate_contains<"CCCCC1NNNN1">(mol);
    //validate_contains<"CCCCC=CC=C1/C(CC=CCCCC(O)=O)CCC1=O">(mol); // FIXME: stereo
    validate_contains<"CCCCCBr">(mol);
    validate_contains<"CCCCCC#C">(mol);
    validate_contains<"CCCCCC(O)C1=CC=CC(OCCC2=NC3=CC=CC=C3C=C2)=C1">(mol);
    //validate_contains<"CCCCCC.CCCCCC.C1CCC2(CC1)CCCCC2">(mol); // FIXME: components
    validate_contains<"CCCCCC1CC(=O)C2=C(C1)CC(C)(C)C3CCC(C)=CC23">(mol);
    validate_contains<"CCCCCC=O">(mol);
    validate_contains<"CCCCCCC1=C2CCNCC2=CC=C1">(mol);
    validate_contains<"CCCCCCC1CCC(=O)O1">(mol);
    validate_contains<"CCCCCCCC1NC2=CC=CC=C2C(=O)C1O">(mol);
    validate_contains<"CCCCCCCCCCCCCCCCCCCCCCN1CCN(CC1)C(=O)C2=CC=C(CC3=NOC(=O)N3)C=C2">(mol);
    validate_contains<"CCCCCCCCCCCCCCCCCCN1CCN(CC2=CC=C(CC3=NOC(=O)N3)C=C2)CC1">(mol);
    validate_contains<"CCCCCCCCCCCCCCCCCCN1CCN(Cc2ccc(CC3=NOC(=O)N3)cc2)CC1">(mol);
    validate_contains<"CCCCCCCCCN1CC(O)C(O)C(O)C1CO">(mol);
    validate_contains<"CCCCCN1C=C(C(=O)C2=CC=CC3=CC=CC=C23)C4=CC=CC=C14">(mol);
    validate_contains<"CCCCCNC(=N)NN=CC1=CNC2=C1C=C(OC)C=C2">(mol);
    validate_contains<"CCCCCOP(O)(O)=O">(mol);
    validate_contains<"CCCCCOS(N)(=O)=O">(mol);
    validate_contains<"CCCCN(CC)C1=NC(C)=NC2=C1C=C(C)N2C3=C(C)C=C(C)C=C3C">(mol);
    validate_contains<"CCCCN1CCCCC1C(=O)NC2=C(C)C=CC=C2C">(mol);
    validate_contains<"CCCCN1CCN(CC1)S(=O)(=O)NCCc1nc([nH]c1-c1ccc(OC)cc1)-c1cccs1">(mol);
    validate_contains<"CCCCO">(mol);
    validate_contains<"CCCCOC1C(=O)N(N(C1=O)C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"CCCN">(mol);
    validate_contains<"CCCN1CCC(CC1)C2=NC3=CC=CC(C(N)=O)=C3N2">(mol);
    validate_contains<"CCCN1CCC(CC1)C2NC3=CC=CC(C(N)=O)=C3N2">(mol);
    validate_contains<"CCCN1CCN(CC1)C2=CC=CC=C2">(mol);
    validate_contains<"CCCNC(=O)C1OC(C(O)C1O)N2C=NC3=C(N)N=CN=C23">(mol);
    validate_contains<"CCCO">(mol);
    validate_contains<"CCCOC1=CC(=CC=C1O)C(C)=O">(mol);
    validate_contains<"CCCOC1=CC=C(C=C1)C2=CC=C(C)C=C2">(mol);
    validate_contains<"CCCOC1=CNC=C1">(mol);
    validate_contains<"CCCOCCC(C(C)C)C1=CC=CC(CCCOC)=C1">(mol);
    validate_contains<"CCCS(=O)(=O)NCC">(mol);
    //validate_contains<"CCN(C(C)C)C(=O)C1=CC(C)=CC(OC[C@H](C)[NH2+]C2=C=CNC=C2)=C1">(mol); // FIXME: stereo
    validate_contains<"CCN(C(C)C)C(=O)c1cc(C)cc(OCC(C)Nc2ccncc2)c1">(mol);
    validate_contains<"CCN(C)C=CCC(P(O)(O)=O)P(O)(O)=O">(mol);
    //validate_contains<"CCN(CC)C(=O)[C@H]1CN(C)C2CC3=CNC4=C3C(=CC=C4)C2=C1">(mol); // FIXME: stereo
    validate_contains<"CCN(CC)CC1=CC(O)=CC(NC2=C3C=CC(Cl)=CC3=NC=C2)=C1">(mol);
    validate_contains<"CCN(CC)CCCC(C)NC1=C2C=C(OC)C=CC2=NC2=C1C=CC(Cl)=C2">(mol);
    validate_contains<"CCN(CC)CCCCNC1=NC=C2C=C(C(N)=NC2=C1)C3=CC(OC)=CC(OC)=C3">(mol);
    validate_contains<"CCN(CC)CCCCNC1=NC=C2C=C(C(N)=NC2=C1)C3=CC=CC=C3">(mol);
}
