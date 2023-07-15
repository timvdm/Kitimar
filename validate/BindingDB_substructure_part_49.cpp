#include "Validate.hpp"

void BindingDB_substructure_part_49(OpenBabel::OBMol &mol)
{
    // SMARTS 2401 - 2450
    validate_contains<"OC(=O)CCC(=O)NC1=CC=C(CC2=CC(O)=C(NC(=O)CCC(O)=O)C=C2)C=C1O">(mol);
    validate_contains<"OC(=O)CCC1=CNC2=C1C=CC=C2">(mol);
    //validate_contains<"OC(=O)CCCC[C@@H]1CCSS1">(mol); // FIXME: stereo
    validate_contains<"OC(=O)CN1C(=O)N(CC2=CC=CC=C2F)C3=C(C4=C(CCCC4)S3)C1=O">(mol);
    validate_contains<"OC(=O)CN1C(=O)NC2=C(C3=C(CCCC3)S2)C1=O">(mol);
    validate_contains<"OC(=O)CNC1CCCCC1">(mol);
    validate_contains<"OC(=O)CSC1=NC(Cl)=CC(NC2=CC=CC=C2)=N1">(mol);
    validate_contains<"OC(=O)N1C=CC(=C1)C2=CC=CC=C2">(mol);
    validate_contains<"OC(=O)P(O)(O)=O">(mol);
    //validate_contains<"OC(=O)[C@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)NC(=O)[C@@H]3C4CCC(CC4)N3C(=O)CCC5CCCC5">(mol); // FIXME: stereo
    validate_contains<"OC(=O)c1ccc(NC(=O)c2ccccc2)cc1">(mol);
    validate_contains<"OC(=O)c1ccccc1O">(mol);
    validate_contains<"OC(CCC=O)CC1CCCC(=O)C1">(mol);
    validate_contains<"OC(NC1=CC=C(OC(F)(F)F)C=C1)C2SCCC2NCC3=C4C=CC=CC4=NC=C3">(mol);
    validate_contains<"OC(NC1=CC=CC=C1)C2=CC3=CC=CC=C3C=C2">(mol);
    //validate_contains<"OC(NN=C/C1=CC=CC=C1)C2=C3C=CC=CC3=NC(=C2)C4=CC=CC=C4">(mol); // FIXME: stereo
    validate_contains<"OC(O)C1=CC=CO1">(mol);
    validate_contains<"OC(O)CC(O)(CC(O)O)C(O)O">(mol);
    validate_contains<"OC(O)CCC(CP(O)(O)=O)C(O)=O">(mol);
    //validate_contains<"OC12CC3CC(C1)CC(C3)(C2)NCC(=O)N4CCC[C@H]4C#N">(mol); // FIXME: stereo
    validate_contains<"OC1=C(C(C2CC2)C3=CC=CC=C3)C(=O)OC=C1">(mol);
    validate_contains<"OC1=C(C2=C(C=CC=C2)C=C1)C3=C(O)C=CC4=C3C=CC=C4">(mol);
    validate_contains<"OC1=C(C=CC=C1C2=CC=CC=C2)C3=CC4=C(N3)C=CC=C4">(mol);
    validate_contains<"OC1=C(CCCOc2ccccc2)C(=O)Oc3ccccc13">(mol);
    validate_contains<"OC1=C(N=C2CCCCCN2C1=O)C(=O)NCc3ccc(F)cc3">(mol);
    validate_contains<"OC1=C(N=C2CCCN2C1=O)C(=O)NCc3ccc(F)cc3">(mol);
    //validate_contains<"OC1=C(N=C2[C@@H](CCCN2C1=O)N3CCCCC3)C(=O)NCc4ccc(F)cc4">(mol); // FIXME: stereo
    //validate_contains<"OC1=C(N=C2[C@H](CCCN2C1=O)N3CCOCC3)C(=O)NCc4ccc(F)cc4">(mol); // FIXME: stereo
    validate_contains<"OC1=C(O)C=C(C=C1)C(=O)CSC2=NCCS2">(mol);
    validate_contains<"OC1=C(O)C=C(C=C1)C(=O)NN=CC2=CC(=CC=C2O)[N+]([O-])=O">(mol);
    validate_contains<"OC1=C(O)C=C(C=CC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"OC1=C2C(=O)CCOC2=CC=C1">(mol);
    validate_contains<"OC1=C2C(=O)OCCC2=CC=C1">(mol);
    validate_contains<"OC1=C2C=CNC2=CC=C1">(mol);
    validate_contains<"OC1=C2CCCC2=C(Cl)C=C1">(mol);
    validate_contains<"OC1=C2N=CNC2=C3NC=CNC(=O)C3=C1">(mol);
    validate_contains<"OC1=C2NC=NC2=C(C=O)C=C1">(mol);
    validate_contains<"OC1=CC(=C)OCC1">(mol);
    validate_contains<"OC1=CC(=CC=C1)N2CCCC2">(mol);
    validate_contains<"OC1=CC(=O)NC2=CC=CC=C12">(mol);
    validate_contains<"OC1=CC(=O)OC=C1">(mol);
    validate_contains<"OC1=CC(O)=C(Cl)C=C1">(mol);
    validate_contains<"OC1=CC2=C(C=C1O)C(=O)C3=CC=CC=C3O2">(mol);
    validate_contains<"OC1=CC2=C(NC=[S]2)C=C1">(mol);
    validate_contains<"OC1=CC2=C(OCCO2)C=C1">(mol);
    validate_contains<"OC1=CC2=CC=CC=C2C(=C1O)C3=C4C=CC=CC4=CC(O)=C3O">(mol);
    validate_contains<"OC1=CC=C(Br)C(O)=C1Br">(mol);
    validate_contains<"OC1=CC=C(C=C1)C(=O)C=CC2=CC=C(Cl)C=C2">(mol);
    validate_contains<"OC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2">(mol);
    //validate_contains<"OC1=CC=C(C=C1)C2=C(Cl)C(=O)OC2=C/C3=CC(Br)=C(O)C(Br)=C3">(mol); // FIXME: stereo
}
