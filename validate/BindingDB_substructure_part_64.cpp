#include "Validate.hpp"

void BindingDB_substructure_part_64(OpenBabel::OBMol &mol)
{
    // SMARTS 3151 - 3200
    validate_contains<"[O-][N+]=O">(mol);
    validate_contains<"[S+2]([O-])([O-])(Cc1[nH0]c(oc1C)c1ccccc1F)CC(=O)NCc1occc1">(mol);
    validate_contains<"c1ccccc1NO">(mol);
    validate_contains<"C1=NC=NC=N1">(mol);
    validate_contains<"C1CN(CCN1)C(C1=CC=CC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"CC(=O)OC1=C(C=CC=C1)C(O)=O">(mol);
    validate_contains<"CN(C)C1=CC=CC=C1">(mol);
    validate_contains<"CNC(=O)C1=CNC2=C1C=CC=C2">(mol);
    validate_contains<"NC(=O)C1=CC=CC=C1">(mol);
    validate_contains<"NC(=S)NN=C">(mol);
    validate_contains<"NC(N)=O">(mol);
    validate_contains<"OC1=C(N=NC2=CC=C(C=C2)S([O-])(=O)=O)C2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"c1ccc2ncccc2c1">(mol);
    //validate_contains<"Br.N1C=CC=C1.C2=CC3=C(C=C2)C=CC=C3">(mol); // FIXME: components
    validate_contains<"BrC1=CC(NC2=NC=NC3=CC=C(NC(=O)C=C)C=C23)=CC=C1">(mol);
    validate_contains<"C1=CC2=C(C=C1)C3=CC=CC=C3C=C2">(mol);
    validate_contains<"C1=CC2=CC=CN=C2C=C1">(mol);
    validate_contains<"C1=CC2=CN=C3C=CC=CC3=C2C=C1">(mol);
    validate_contains<"C1=CC=C(C=C1)N2C=CC=N2">(mol);
    validate_contains<"C1CC2=CC=CC=C2CN1">(mol);
    validate_contains<"C1CCOC1">(mol);
    validate_contains<"C1CN(CCN1)C2=CC=CC=C2">(mol);
    validate_contains<"C1COC1">(mol);
    validate_contains<"C=S">(mol);
    //validate_contains<"CC(=C/C(=O)NC1=C(C=CC=C1)C(O)=O)C2=CC3=CC=CC=C3C=C2">(mol); // FIXME: stereo
    validate_contains<"CC(=O)OC1=CC=CC=C1C(O)=O">(mol);
    validate_contains<"CC(C)C(NC(=O)N(C)CC1=CSC(=N1)C(C)C)C(=O)NC(CC(O)C(CC2=CC=CC=C2)NC(=O)OCC3=CN=CS3)CC4=CC=CC=C4">(mol);
    //validate_contains<"CC(C)[C@H](CO)Nc1nc(Nc2ccc(C(O)=O)c(Cl)c2)c3ncn(C(C)C)c3n1">(mol); // FIXME: stereo
    //validate_contains<"CC(C)[C@H](NC(=O)OCc1ccccc1)C(=O)N[C@@H](Cc2ccccc2)C(O)[C@H](Cc3ccccc3)NC(=O)[C@@H](NC(=O)OCc4ccccc4)C(C)C">(mol); // FIXME: stereo
    validate_contains<"CC(N(C)C)c1cccc(O)c1">(mol);
    //validate_contains<"CC.CCC">(mol); // FIXME: components
    validate_contains<"CC1(C)CC(=O)C2=C(O1)C=CC=C2">(mol);
    validate_contains<"CC12CCC3C(CCC4=CC(O)=CC=C34)C1CCCC2O">(mol);
    validate_contains<"CC1=CC2=C(C=CC=C2)C=N1">(mol);
    validate_contains<"CC1=CC=CC=C1C">(mol);
    validate_contains<"CC1CCN(C)CC1">(mol);
    validate_contains<"CCC(C1=CC=CC=C1)C2=C(O)C3=C(OC2=O)C=C(OC)C=C3">(mol);
    //validate_contains<"CCCNC(=O)[C@@H]1CCCN1C(=O)[C@H]([NH3+])C(C2=CC=CC=C2)C3=CC=CC=C3">(mol); // FIXME: stereo
    validate_contains<"CCCOC1=CC=C(F)C=C1">(mol);
    //validate_contains<"CCC[C@@]2(CCC1=CC=CC=C1)CC(O)=C([C@H](CC)C3=CC(NS(=O)(=O)C4=NC=C(C=C4)C(F)(F)F)=CC=C3)C(=O)O2">(mol); // FIXME: stereo
    validate_contains<"CCN1C2=C(C=CC=N2)C(=O)N(C)C3=C1N=C(C=C3)N4CCOCC4">(mol);
    validate_contains<"CCNC(=O)c1noc(-c2cc(C(C)C)c(O)cc2O)c1-c3ccc(CN4CCOCC4)cc3">(mol);
    validate_contains<"CCOC(=O)CN1CCN(CC1)c2ccc(cc2NC(=O)c3cccc4ccccc34)-c5ccccc5">(mol);
    validate_contains<"CN(C)CC1=Nc2c(sc3ccccc23)C(=O)N1">(mol);
    validate_contains<"CN1CCC(CC1)OC(C2=CC=CC=C2)C3=CC=CC=C3">(mol);
    validate_contains<"CN1CCCC1CC2=CC=CC=C2">(mol);
    validate_contains<"CNC(=O)C1=CC=CC=C1">(mol);
    validate_contains<"CNC1=SC=CN1">(mol);
    validate_contains<"COC1=CC(OC)=C2C(NC3=C4OCOC4=CC=C3Cl)=NC=NC2=C1">(mol);
    //validate_contains<"C[C@@H]1C[C@H]2O[C@@H]2C=CC=CC(=O)Cc3c(Cl)c(O)cc(O)c3C(=O)O1">(mol); // FIXME: stereo
}
