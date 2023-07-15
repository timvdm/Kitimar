#include "Validate.hpp"

void BindingDB_substructure_part_12(OpenBabel::OBMol &mol)
{
    // SMARTS 551 - 600
    validate_contains<"*C1=CC2=C(N1)N=CC(*)=C2">(mol);
    validate_contains<"*N1C=CN=C1">(mol);
    validate_contains<"*NC(*)C(=O)N*">(mol);
    //validate_contains<"*NC(*)C(=O)N*.NC(=N)C1=CC=CC=C1">(mol); // FIXME: components
    validate_contains<"*NC1=CC(N*)=NC=N1">(mol);
    //validate_contains<"Br.C1=CC2=C(C=C1)C=CC=C2">(mol); // FIXME: components
    //validate_contains<"Br.N1C=CC=C1">(mol); // FIXME: components
    validate_contains<"BrC(C1CCCCC1)C2=CNC=C2">(mol);
    validate_contains<"BrC1=C(Br)C=CC=C1">(mol);
    validate_contains<"BrC1=CC(C=C2OCNC2=O)=CC(Br)=C1">(mol);
    validate_contains<"BrC1=CC(NC2=NC=NC3=C2C=C(NC(=O)C=C)C=C3)=CC=C1">(mol);
    validate_contains<"BrC1=CC(NC2=NC=NC3=CC(NC(=O)C=C)=CC=C23)=CC=C1">(mol);
    validate_contains<"BrC1=CC=CC(NC2=NC=NC3=CC(NC(=O)C=C)=CC=C23)=C1">(mol);
    validate_contains<"BrC1=CC=CN1">(mol);
    validate_contains<"BrC1=NN=C(SCC(=O)NC2=CC=CC=C2)N1C3=CC=CC=C3">(mol);
    validate_contains<"BrC1CCCC1">(mol);
    validate_contains<"BrC=CC1=CNC(=O)NC1=O">(mol);
    //validate_contains<"BrCC1=CC(Br)=CC=C1.CC=C/C=C2CC=CC=C2">(mol); // FIXME: components
    validate_contains<"BrCC1CCCCC1">(mol);
    validate_contains<"Brc1cc(ccc1)Nc2ncnc3c2cc(c(c3)NC)N">(mol);
    //validate_contains<"Brc1ccc(NC(=N)N[C@@H]2C[C@@H]2c3ccccc3)nc1">(mol); // FIXME: stereo
    validate_contains<"Brc1cccc(Nc2ncnc3ccccc23)c1">(mol);
    validate_contains<"C#CC1=CC(=CC(=C1)C#C)C#C">(mol);
    validate_contains<"C#N">(mol);
    validate_contains<"C(C)(=O)N(C(C(=O)N(C(C)(COC(=O)OC)C)[H])C(C)C)C">(mol);
    validate_contains<"C(C1CCCC1)C2=CN=C3C2=CC4=C5C(C=CC=C35)=CC=C4">(mol);
    validate_contains<"C(C1CCCC2=C1C=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"C(C1CNC2=CC=CC=C12)C3=CCNC3">(mol);
    validate_contains<"C(CC1=CC=C2C=CC=CC2=C1)C3=CC4=CC=CC=C4C=C3">(mol);
    validate_contains<"C(CC1=CC=CC2=CC=CC=C12)N3C=CC=C3">(mol);
    validate_contains<"C(CC1=CC=CC=C1)N2CCCCC2">(mol);
    validate_contains<"C(CC1=CC=CN=C1)C2=CC=CC=C2">(mol);
    validate_contains<"C(CC1=CNC2=CC=CC=C12)N3C=CC=C3">(mol);
    validate_contains<"C(CC1CC1)CC2CC2">(mol);
    validate_contains<"C(CC1CCCC1)NCCN2CCN(C2)C3=CC=CC=C3">(mol);
    validate_contains<"C(CN(CCC1=CC=CC=C1)CC2=CC=CC=C2)CC3=CC=CC=C3">(mol);
    validate_contains<"C(N1CC2=CC=CC=C2C1)C3=CC=CC=C3">(mol);
    validate_contains<"C(OC1=CC=C(C=C1)C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"C(Oc1ccccc1)C2CCCN2Sc3ccc4N(Cc5ccccc5)CCc4c3">(mol);
    validate_contains<"C(c1ccncc1)c2nnc(Nc3ccc(cc3)-c4ccccc4)c5ccccc25">(mol);
    //validate_contains<"C.C.CC.CC.C1CCCC1.CCCCCCCC.CCC(C)CC2(C)CCCC2">(mol); // FIXME: components
    //validate_contains<"C.CC.OCCCNc1cncc(c1)-c2cncc(Nc3cccc(Cl)c3)n2">(mol); // FIXME: components
    //validate_contains<"C.CCC(CC)C(=O)C1CCC2CC3CNC4=C3C(=CC=C4)C2=C1">(mol); // FIXME: components
    //validate_contains<"C.N1C=CC2=C1C=C3C=CC4=CC=CC=C4C3=C2">(mol); // FIXME: components
    validate_contains<"C1(=CC(=NC(=C1)CC)CO)CCCC">(mol);
    //validate_contains<"C1(=CC(=NC(=C1)CC)CO)CCCC.[H]">(mol); // FIXME: components
    validate_contains<"C1*CC2=C1C3=C(*C4=CC=CC=C34)C5=C2C6=C(*5)C=CC=C6">(mol);
    validate_contains<"C1=C(N2N=CC=CC2=N1)C3=NC=NC=C3">(mol);
    validate_contains<"C1=CC2=C(C=C1)C3=C(C=C2)C4=CN=CC=C4C=C3">(mol);
    validate_contains<"C1=CC2=C(C=C1)C3=C(C=CC=C3)C=C2">(mol);
}
