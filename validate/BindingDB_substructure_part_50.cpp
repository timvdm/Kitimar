#include "Validate.hpp"

void BindingDB_substructure_part_50(OpenBabel::OBMol &mol)
{
    // SMARTS 2451 - 2500
    validate_contains<"OC1=CC=C(C=C1)C2=CC(=NO2)C(=O)NN=C">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2=CC(=NO2)C(=O)NN=CC3=CC=CC=C3">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2=COC3=CC(O)=CC(O)=C3C2=O">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2=NC(NCC3=CC=CC=C3)=CC(=C2)C4=CC=CC=C4">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2=NC(NCCC3=CC=CC=C3)=CC(=C2)C4=CC=CC=C4">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2CC(=O)C3=C(C2)C=C(O)C=C3O">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2CC(=O)C3=C(O2)C=C(O)C=C3O">(mol);
    validate_contains<"OC1=CC=C(C=C1)C2OC3=CC=C(O)C=C3C4CCCC24">(mol);
    validate_contains<"OC1=CC=C(C=C1)C=CC(=O)C2=C(O)C=C(O)C=C2O">(mol);
    validate_contains<"OC1=CC=C(C=C1)C=CC2=CC(O)=CC(O)=C2">(mol);
    validate_contains<"OC1=CC=C(C=C2SC(=N)NC2=O)C=C1O">(mol);
    validate_contains<"OC1=CC=C(C=NNC(=O)C2=NOC(=C2)C3=CC=C(O)C=C3)C=C1">(mol);
    validate_contains<"OC1=CC=C(CC=C)C=C1C2=C(O)C=CC(CC=C)=C2">(mol);
    validate_contains<"OC1=CC=C(CCC2=CC=CC=C2)C=C1">(mol);
    validate_contains<"OC1=CC=C(CN2C3=CC=CC=C3C4=C2CN(CC4)C(=O)C5=CC=CC=C5)C=C1">(mol);
    validate_contains<"OC1=CC=C(Cl)C(O)=C1">(mol);
    validate_contains<"OC1=CC=C(Cl)C2=CC=CN=C12">(mol);
    validate_contains<"OC1=CC=CC(O)=C1">(mol);
    validate_contains<"OC1=CC=CC2=C1C(O)=CC=C2">(mol);
    validate_contains<"OC1=CC=CC=C1C([O-])=O">(mol);
    validate_contains<"OC1=CC=CC=C1C2=NC3=C(N2)C=CC=C3">(mol);
    validate_contains<"OC1=CC=NC2=C1C=C(F)C=C2">(mol);
    validate_contains<"OC1=CN=C(Cl)C2=C1C=CC=C2">(mol);
    validate_contains<"OC1=CN=CC2=C1C=CC=C2">(mol);
    validate_contains<"OC1=Cc2cc(O)c(O)cc2OC1=O">(mol);
    validate_contains<"OC1=NC2=C(C=CC=C2)C(O)=N1">(mol);
    validate_contains<"OC1=NC=CC2=C1C=CC=C2">(mol);
    validate_contains<"OC1=NC=CC=C1">(mol);
    validate_contains<"OC1=NOC=C1C(=O)NN=CC2=CC=CC=C2">(mol);
    validate_contains<"OC1C(OC2=C(C(O)=CC=C2)C1=O)C3=CC(O)=C(O)C=C3">(mol);
    validate_contains<"OC1CCC(COP(O)(O)=O)C1O">(mol);
    validate_contains<"OC1CCC(N1)C(O)=O">(mol);
    validate_contains<"OC1CCCCC1N2CCC(CC2)C3=CC=CC=C3">(mol);
    validate_contains<"OC1CCCN(C1)CC2=CC=CC=C2">(mol);
    validate_contains<"OC1CCCO1">(mol);
    validate_contains<"OC1CCOC(=O)C1">(mol);
    validate_contains<"OC1CCOC1O">(mol);
    validate_contains<"OC1COCC1O">(mol);
    validate_contains<"OC1NC(O)C2=CC=CC=C2N1">(mol);
    validate_contains<"OC=CC1=NC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"OC=O">(mol);
    //validate_contains<"OC=O.CCCC1=NN(C=C1)C2=CC=CC=C2">(mol); // FIXME: components
    validate_contains<"OCC(=O)CC1=CC=CC2=C1C=CC=C2">(mol);
    validate_contains<"OCC(=O)COP(O)(O)=O">(mol);
    validate_contains<"OCC(NC(=O)C(Cl)Cl)C(O)C1=CC=C(C=C1)[N+]([O-])=O">(mol);
    validate_contains<"OCC(O)=O">(mol);
    validate_contains<"OCC(O)C(O)C1=CNC2=CC=CC=C12">(mol);
    validate_contains<"OCC(O)CO[Na]">(mol);
    validate_contains<"OCC1=CC=C(O1)C=O">(mol);
    validate_contains<"OCC1CC(C=C1)N2C=CC(NO)=NC2=O">(mol);
}
