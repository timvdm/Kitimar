#include "Validate.hpp"

void BindingDB_substructure_part_61(OpenBabel::OBMol &mol)
{
    // SMARTS 3001 - 3050
    validate_contains<"FC1=C(Br)C=CC=C1Br">(mol);
    validate_contains<"FC1=CC=C(CN2C(=O)C(=O)C3=CC(=CC=C23)S(=O)(=O)N4CCCC4COC5=C(F)C=C(F)C=C5)C=C1">(mol);
    validate_contains<"FC1=CC=C(SC2=NN3C=NC(=O)C(=C3C=C2)C4=C(Cl)C=CC=C4Cl)C(F)=C1">(mol);
    validate_contains<"Fc1cc(cc(c1)C(=O)Nc2ccc(Cl)c(COc3cccnc3)c2)N4CCOCC4">(mol);
    //validate_contains<"Fc1ccc(cc1)[C@@H]2CCNC[C@H]2COc3ccc4OCOc4c3">(mol); // FIXME: stereo
    validate_contains<"N#CC(C#N)=C1C=CC=CC=C1">(mol);
    validate_contains<"N#Cc1ccc2cc(ccc2c1)-c3cccnc3">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=C3C=CC=CC3=NC=N2">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"N(C1=CC=CC=C1)C2=NC=NC3=CC=CC=C23">(mol);
    validate_contains<"N1C2=C(SC3=C1C=CC=C3)C=CC=C2">(mol);
    validate_contains<"N1C2=CC=CC=C2SC3=C1C=CC=C3">(mol);
    validate_contains<"N1C=C2C=CC=CC2=N1">(mol);
    validate_contains<"N1C=CC2=C1C=CC=N2">(mol);
    validate_contains<"N1C=CC2=C1N=CN=C2">(mol);
    validate_contains<"N1C=CC2=C1NC=C2">(mol);
    validate_contains<"N1C=CC2=CC3=C(C=CC=C3)C=C12">(mol);
    validate_contains<"N1C=NC2=C1N=CC=C2">(mol);
    validate_contains<"N1N=NC2=C1C=CC=C2">(mol);
    validate_contains<"NC(=N)C1=CC2=C(NC(=C2)C3=CC=CC(C4=CC=CC=C4)=C3O)C=C1">(mol);
    validate_contains<"NC(=N)C1=CC=C(C=C1)C2=CC=CO2">(mol);
    validate_contains<"NC(=N)c1ccc(NC(=O)c2ccccc2O)cc1">(mol);
    validate_contains<"NC(=O)C1=C(N)C=CS1">(mol);
    validate_contains<"NC(=O)C1=C(S)C2=C(N1)C=CC=C2">(mol);
    validate_contains<"NC(=O)C1=CN=C(N)S1">(mol);
    validate_contains<"NC(=O)C1=CNC=N1">(mol);
    validate_contains<"NC(=O)C1=COC2=CC=CC=C12">(mol);
    validate_contains<"NC(=O)CN1C=CC=C1C(N)=O">(mol);
    validate_contains<"NC(=O)NC1=CC=CC=C1">(mol);
    validate_contains<"NC(=O)NN=C">(mol);
    validate_contains<"NC(=S)NN=Cc1ccc(O)cc1O">(mol);
    validate_contains<"NC(=[NH2+])c1ccc2[nH]c(nc2c1)C(=O)N(O)Cc3cccc(Oc4ccccc4)c3">(mol);
    validate_contains<"NC1=C2CCCCC2=NC3=CC=CC=C13">(mol);
    validate_contains<"NC1=C2NCN(C3OC(O)C(O)C3O)C2=NC=C1">(mol);
    validate_contains<"NC1=CC(=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"NC1=CC(C2=CC=C(C=C2)S(O)(=O)=O)=C(O)C(=C1)C(O)=O">(mol);
    validate_contains<"NC1=CC(Cl)=NC(N)=N1">(mol);
    validate_contains<"NC1=CC(N)=NC(N)=N1">(mol);
    validate_contains<"NC1=CC2=CC=C(N)N=C2N=C1">(mol);
    validate_contains<"NC1=CC=C(N)C=C1">(mol);
    validate_contains<"NC1=CC=CC=C1N">(mol);
    //validate_contains<"NC1=NC(=O)C2=C(NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)=N2)N1">(mol); // FIXME: stereo
    //validate_contains<"NC1=NC2=CC=C(NC(=O)[C@@H]3CCCN3C(=O)[C@@H](CC4CCCCC4)NS(=O)(=O)CC5=CC=CC=C5)C=C2S1">(mol); // FIXME: stereo
    validate_contains<"NC1=NC=NC2=C1C=NN2">(mol);
    validate_contains<"NC1=NC=NC2=C1N=CN2CC3=CC=CC=C3">(mol);
    validate_contains<"NC1=[S]C2=C(N1)C=CC=C2">(mol);
    validate_contains<"NCC(=O)NCC(O)=O">(mol);
    validate_contains<"NCC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2N">(mol);
    validate_contains<"NCC1=CC=C(C=C1)C2=CC=CC=C2">(mol);
    validate_contains<"NCC1=CC=C(NC(=O)C2=CC=C3C=C(C=CC3=C2)C(N)=N)C=C1">(mol);
}
