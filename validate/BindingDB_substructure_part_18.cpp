#include "Validate.hpp"

void BindingDB_substructure_part_18(OpenBabel::OBMol &mol)
{
    // SMARTS 851 - 900
    validate_contains<"CC(=O)OCC1OC(C(OC(=S)NCC2=CC=CC=C2)C(OC(C)=O)C1OC(C)=O)N3CCCCC3">(mol);
    validate_contains<"CC(=O)OCCBr">(mol);
    validate_contains<"CC(=O)OCC[N+](C)(C)C">(mol);
    validate_contains<"CC(C(=O)NC1=CC=C(I)C=C1C(N)=O)=C(C)C2=CC=CC=C2">(mol);
    validate_contains<"CC(C(N)C(C(O)=O)C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"CC(C(O)=O)C(=O)NC1(C)CCCC1">(mol);
    validate_contains<"CC(C)(C#N)C1=CC(=CC(CN2C=NC=N2)=C1)C(C)(C)C#N">(mol);
    validate_contains<"CC(C)(C)C">(mol);
    validate_contains<"CC(C)(C)C(=O)N1CC(=CC1C2=CC=CC=C2)C3=CC(F)=CC=C3F">(mol);
    validate_contains<"CC(C)(C)C1=CC=C(C=C1)C(O)=O">(mol);
    validate_contains<"CC(C)(C)C1=CC=C(C=C1)C2=NN=C(O2)C3=C(OC4=CN=CC(Cl)=C4)C=CC(=C3)N(=O)=O">(mol);
    validate_contains<"CC(C)(C)C1=CN=C(CSC2=CN=C(NC(=O)C3CCNCC3)S2)O1">(mol);
    validate_contains<"CC(C)(C)C1=NC2=C(C=CC=C2C(=C1)N=NCC3=CC=CC=C3)C(C)(C)C">(mol);
    validate_contains<"CC(C)(C)C1=NC2=C3C=CC(F)=CC3=C4C(=O)NC=CC4=C2[N]1">(mol);
    validate_contains<"CC(C)(C)C1CCCCC1">(mol);
    validate_contains<"CC(C)(C)N">(mol);
    validate_contains<"CC(C)(C)N1CC(=O)N=CC2=CC(N)=CC=C12">(mol);
    validate_contains<"CC(C)(C)N1N=C(C2=C1N=CN=C2N)C3=CC=C(Cl)C=C3">(mol);
    validate_contains<"CC(C)(C)N1N=C(C2=CC=C(Cl)C=C2)C3=C1N=CN=C3N">(mol);
    validate_contains<"CC(C)(C)N1N=C(C2=CC=CC3=C2C=CC=C3)C4=C1N=CN=C4N">(mol);
    validate_contains<"CC(C)(C)N1NC(C2=CC=CC3=C2C=CC=C3)C4=C(N)NCNC14">(mol);
    validate_contains<"CC(C)(C)NC(=O)C1CC2CCCCC2CN1CC(O)C(CC3=CC=CC=C3)NC(=O)C(CC(N)=O)NC(=O)C4=CC=C5C=CC=CC5=N4">(mol);
    validate_contains<"CC(C)(C)NC(=O)C1CC2CCCCC2CN1CC(O)C(CC3=CC=CC=C3)NC(=O)C(CC(N)=O)NC(=O)C4=NC5=CC=CC=C5C=C4">(mol);
    validate_contains<"CC(C)(C)NC(=O)CNC=O">(mol);
    //validate_contains<"CC(C)(C)NC(=O)[C@@H]1CC2CCCCC2CN1C[C@@H](O)[C@H](CC3=CC=CC=C3)NC(=O)[C@H](CC(N)=O)NC(=O)C4=CC=C5C=CC=CC5=N4">(mol); // FIXME: stereo
    //validate_contains<"CC(C)(C)NC(=O)[C@@H]1CN(CC2=CC=CN=C2)CCN1C[C@@H](O)C[C@@H](CC3=CC=CC=C3)C(=O)N[C@@H]4[C@H](O)CC5=CC=CC=C45">(mol); // FIXME: stereo
    validate_contains<"CC(C)(C)NC(CC1=CC=C2OCOC2=C1)C(C)(C)C">(mol);
    validate_contains<"CC(C)(C)OC(=O)NC1(CCCC1)C2=NC3=NC=NC=C3N=C2">(mol);
    validate_contains<"CC(C)(C)c1cc(NC(=O)Nc2cccc(Cl)c2Cl)n(n1)-c3ccccc3">(mol);
    //validate_contains<"CC(C)(C)c1ccc(cc1)S(=O)(=O)N[C@@H](Cc2cccc(c2)[CH+](N)N)C(=O)N3CCN(CC3)S(C)(=O)=O">(mol); // FIXME: stereo
    validate_contains<"CC(C)(CCl)CC1=CC(=C(N)C(=C1)C(C)(C)Cl)C(C)(C)N">(mol);
    validate_contains<"CC(C)(COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1OP(O)(O)=O)N2C=NC3=C2N=CN=C3N)C(O)C(=O)NCCC(=O)NCCSCC(=O)NCCC4=CNC5=CC=CC=C45">(mol);
    validate_contains<"CC(C)(N)C1=CC=CC=C1">(mol);
    validate_contains<"CC(C)(N)CC1=CC=CC=C1">(mol);
    validate_contains<"CC(C)=CCC(O)C1=CC(=O)C2=C(O)C=CC(O)=C2C1=O">(mol);
    validate_contains<"CC(C)=CCc1cc2c(cc1O)oc1cc(O)c(O)c(CC=C(C)C)c1c2=O">(mol);
    //validate_contains<"CC(C)=N/O">(mol); // FIXME: stereo
    validate_contains<"CC(C)=S">(mol);
    validate_contains<"CC(C)C(=O)N(C)C1CCC2(C)C(CCC3C2CCC4(C)CCC(C(C)N(C)C)C34)C1">(mol);
    validate_contains<"CC(C)C(=O)OC(C)(C)C">(mol);
    validate_contains<"CC(C)C(N)C(=O)N1CCCC1">(mol);
    validate_contains<"CC(C)C(NC(=O)N(C)CC1=CNC(=[S]1)C(C)C)C(=O)NC(CC(O)C(CC2=CC=CC=C2)NC(=O)OCC3=CNC=[S]3)CC4=CC=CC=C4">(mol);
    validate_contains<"CC(C)C(NC(C)=O)C(=O)NCCC(O)=O">(mol);
    validate_contains<"CC(C)C(O)C1(NC(=O)C(C)C1O)C(=O)SCC(NC(C)=O)C(O)=O">(mol);
    //validate_contains<"CC(C)C.ClC1=CC=CC=C1">(mol); // FIXME: components
    validate_contains<"CC(C)C1=C(SC2=CC(Cl)=CC(Cl)=C2)N(CC3=CC=NC=C3)C(COC(N)=O)=N1">(mol);
    validate_contains<"CC(C)C1=CC(CC2C(=O)N(CNC(=O)C3=CC=CC=C3)C4=[S]C=C(OCC5=CC=CC(=C5)C(O)=O)C=C24)=CC(C(C)C)=C1C">(mol);
    validate_contains<"CC(C)C1=CC(CCC[N+](O)=O)=CC(C(C)C)=C1F">(mol);
    validate_contains<"CC(C)C1=CC=C(NC(C)=O)C=C1">(mol);
    validate_contains<"CC(C)C1=CC=C2CC2=C1">(mol);
}
