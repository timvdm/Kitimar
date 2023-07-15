#include "Validate.hpp"

void BindingDB_substructure_part_34(OpenBabel::OBMol &mol)
{
    // SMARTS 1651 - 1700
    //validate_contains<"CO.CC12CCCC1C3CCC4=CC(O)=CC=C4C3CC2">(mol); // FIXME: components
    validate_contains<"COC(=O)C(=O)NC(C)(C)C1=NC(C(=O)NCc2ccc(F)cc2)=C(O)C(=O)N1C">(mol);
    validate_contains<"COC(=O)C(CC1=CC=CC=C1)NS(=O)(=O)C2=CC=CC=C2">(mol);
    validate_contains<"COC(=O)C(N1CCC2=C(C1)C=CS2)C3=CC=CC=C3Cl">(mol);
    validate_contains<"COC(=O)C1=C(C=CC=C1)C(=O)OC">(mol);
    validate_contains<"COC(=O)C1=C(C=CC=C1)C(O)=O">(mol);
    validate_contains<"COC(=O)C1=CC(=CC=C1)C(=C)C2=CC(=CC=C2)C(=O)OC">(mol);
    validate_contains<"COC(=O)C1=CC(CO)=CS1">(mol);
    validate_contains<"COC(=O)C1=CC=C(N1C)S(=O)(=O)NC2=CC=C3COC(=O)C3=C2">(mol);
    validate_contains<"COC(=O)C1=CC=CC(Cl)=C1C">(mol);
    validate_contains<"COC(=O)C1=CC=CC=C1C(=O)OC">(mol);
    validate_contains<"COC(=O)C1CCCN1C(=O)CNCC2=CC=CC=C2">(mol);
    validate_contains<"COC(=O)C=C1NC2=CC=CC=C2NC1=O">(mol);
    validate_contains<"COC(=O)C=C1NC2=CC=CC=C2OC1=O">(mol);
    validate_contains<"COC(=O)CC(C1=CC=CC=C1)P(O)(O)=O">(mol);
    validate_contains<"COC(=O)CCCCC1CCSS1">(mol);
    validate_contains<"COC(=O)CNC(=O)C1=CC(=O)OC2=CC=CC=C12">(mol);
    validate_contains<"COC(=O)CNC(=O)C1CCCN1">(mol);
    validate_contains<"COC(=O)NC1=CC2=C(C=C1)C3=CC=C(C=C3O2)S(=O)(=O)NC(C(C)C)C(O)=O">(mol);
    validate_contains<"COC(=O)NC1=CC2=C(OC3=CC(=CC=C23)S(=O)(=O)NC(C(C)C)C(O)=O)C=C1">(mol);
    validate_contains<"COC(=O)OCC1=CC=C(OP(O)(O)=O)C=C1">(mol);
    //validate_contains<"COC(CO)CC(O)CO.OCC1CC(O)C(O)CO1.OCC2CC(O)C(O)CO2.OCC3CC(O)C(O)CO3.OCC4CC(O)C(O)CO4.OCC5CC(O)C(O)CO5">(mol); // FIXME: components
    //validate_contains<"COC1(C)C(O)CC2=C(C=O)C=CC(C=C(/C)C)=C12">(mol); // FIXME: stereo
    //validate_contains<"COC12CC[C@@]3(CC1C(C)(C)O)[C@H]4Cc5ccc(O)c6OC2[C@]3(CCN4CC7CC7)c56">(mol); // FIXME: stereo
    validate_contains<"COC1=C(C(O)C2=C(N)N=C(NC3CCN(CC3)S(C)(=O)=O)N=C2)C(F)=C(F)C=C1">(mol);
    validate_contains<"COC1=C(C=C2C(=C1)CC(C2=O)CC3=CC=NC=C3)OC">(mol);
    validate_contains<"COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC">(mol);
    validate_contains<"COC1=C(C=CC=C1)N2CCN(CCCCN3C(=O)CCC3=O)CC2">(mol);
    validate_contains<"COC1=C(F)C(F)=CC(F)=C1F">(mol);
    validate_contains<"COC1=C(O)C=C2CCC(C)(C)OC2=C1">(mol);
    //validate_contains<"COC1=C(O)C=CC(C=C2/CCC(CNC3=CC=CC=C3)C2=O)=C1">(mol); // FIXME: stereo
    //validate_contains<"COC1=C(OC)C(=O)C(C)C(CC=C(/C)CCC=C(/C)CCC=C(C)C)C1O">(mol); // FIXME: stereo
    validate_contains<"COC1=C(OC)C(OC)=C2C=CC(=O)OC2=C1">(mol);
    validate_contains<"COC1=C(OC)C2=C(CC3N(CCC4=C3C=C5OCOC5=C4)C2)C=C1">(mol);
    validate_contains<"COC1=C(OC)C=C(C=C1)C1=NNC(=O)CC1">(mol);
    validate_contains<"COC1=C(OC)C=C(C=CC(=O)NC2=CC=CC=C2C(O)=O)C=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(=C1)N=CN=C2N1CCN(CC1)C(O)=O">(mol);
    validate_contains<"COC1=C(OC)C=C2C(=O)C(CC3CCNCC3)CC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(NC3=CC(C)=CC=C3)=NN=CC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(NC3=CC(Cl)=C(F)C=C3)=NC=NC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(NC3=CC=CC=C3)=C(C=NC2=C1)C#N">(mol);
    validate_contains<"COC1=C(OC)C=C2C(NC3=CC=CC=C3)=CC=NC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2C(OC3=CC=CC=C3)=CC=NC2=C1">(mol);
    validate_contains<"COC1=C(OC)C=C2N=C3C(N)=NC4=C(C)C(N)=CC=C4C3=CC2=C1">(mol);
    validate_contains<"COC1=C(OCC2CCN(C)CC2)C=C3N=CN=C(NC4=CC=CC=C4)C3=C1">(mol);
    validate_contains<"COC1=C(OCCCN2CCOCC2)C=C3C(NC4=CC=CC(=C4)C#C)=NC=NC3=C1">(mol);
    validate_contains<"COC1=C2C=CC3=C(NC(=O)NC3=O)C2=CC=C1">(mol);
    validate_contains<"COC1=C2C=CC=CC2=CC3=CC=CC=C13">(mol);
    validate_contains<"COC1=CC(=CC(OC)=C1OC)C2=CN(C)C=N2">(mol);
    validate_contains<"COC1=CC(Br)=CC(C(=O)NC2=CC=C(F)C=C2)=C1OC">(mol);
}
