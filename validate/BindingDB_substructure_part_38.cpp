#include "Validate.hpp"

void BindingDB_substructure_part_38(OpenBabel::OBMol &mol)
{
    // SMARTS 1851 - 1900
    validate_contains<"ClC1CC(NS(=O)(=O)C2=CC=C(Cl)C=C2)=CC(Cl)=C1">(mol);
    validate_contains<"ClC1CCCCC1Br">(mol);
    //validate_contains<"ClC1CCCCC1Br.ClC2CCCC[C@@H]2Br.ClC3CCCC[C@@H]3Br">(mol); // FIXME: components
    //validate_contains<"ClC1CCCC[C@@H]1Br">(mol); // FIXME: stereo
    //validate_contains<"ClCC(=O)N=C1/NC2CS(=O)(=O)CC2S1">(mol); // FIXME: stereo
    validate_contains<"ClN1C=CC(CC2=CC=CC(Cl)=C2)=C1">(mol);
    validate_contains<"Cl[Ti]Cl">(mol);
    validate_contains<"Clc1ccc(NC(=O)c2scnc2CCc3cnoc3)cc1">(mol);
    validate_contains<"Clc1ccc(NC(=O)c2snnc2c3ccccc3)cc1">(mol);
    validate_contains<"Clc1ccc(cc1)C(=O)N2CCN(CC2)c3ccc(c4ncccc34)N(=O)=O">(mol);
    validate_contains<"Cn1cc(C2=C(C(=O)NC2=O)c3cn(CCCNC(N)=N)c4ccccc34)c5ccccc15">(mol);
    validate_contains<"FC(F)(F)C1=CC(=CC=C1)C2=CC3=C(NC=N3)C=C2">(mol);
    validate_contains<"FC(F)(F)C1=CC(=CC=C1)N2CCN(CCN3C(=O)NC4=C3C=CC=C4)CC2">(mol);
    validate_contains<"FC(F)(F)C1=CC(=CC=C1)N2CCNCC2">(mol);
    validate_contains<"FC(F)(F)C1=CC(NC(=O)C2=C(NCC3=CC=NC=C3)C=CC=C2)=CC=C1">(mol);
    validate_contains<"FC(F)(F)C1=CC(NC(=O)C2=CC(CC3=C(C=CC=N3)C4=NC=NC=C4)=CC=C2)=CC=C1">(mol);
    validate_contains<"FC(F)(F)C1=CC(NC(=O)C2=CC=CC=C2NCC3=CC=NC=C3)=CC=C1">(mol);
    validate_contains<"FC(F)(F)C1=CC(NC2=CC=NC(NC3=CC(NC(=O)C4CC4)=CC=C3)=N2)=CC=C1">(mol);
    validate_contains<"FC(F)(F)C1=CC2=C(C=CC=C2C(CN=NC3=CC=CC=C3)=C1)C(F)(F)F">(mol);
    validate_contains<"FC(F)(F)C1=CC=C2SC3=C(NC2=C1)C=CC=C3">(mol);
    validate_contains<"FC(F)(F)C1=NC2=C(C=CC=C2C(F)(F)F)C=C1">(mol);
    validate_contains<"FC(F)(F)OC1=CC=C(NC(=O)C2=C(NCC3=C4C=CC=CC4=NC=C3)C=CS2)C=C1">(mol);
    validate_contains<"FC(F)(F)OC1=CC=C(NC(=O)C2=C(SCC3=C4C=CC=CC4=NC=C3)C=CS2)C=C1">(mol);
    validate_contains<"FC(F)(F)c1ccc(cc1)c2nc(Nc3ccc(Cl)c(Cl)c3)nc(n2)c4ccc(cc4)C(F)(F)F">(mol);
    validate_contains<"FC1=C(Cl)C=C(NC2=C(C=NC3=CC=C(NC(=O)C=C)C=C23)C#N)C=C1">(mol);
    validate_contains<"FC1=C(Cl)C=C(NC2=NC=NC3=CC(OCCN4CCOCC4)=C(NC(=O)C=C)C=C23)C=C1">(mol);
    validate_contains<"FC1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"FC1=C(Cl)N=C(Cl)C(=C1)C#N">(mol);
    validate_contains<"FC1=CC(Br)=CC(F)=C1OC(Cl)=S">(mol);
    validate_contains<"FC1=CC=C(C=C1)C(COC=O)C2CCC2">(mol);
    validate_contains<"FC1=CC=C(C=C1)C1=C(NC(=N1)C1=CC=CC=C1)C1=CC=NC=C1">(mol);
    validate_contains<"FC1=CC=C(C=C1)C2=C(NC=N2)C3=CC=NC=C3">(mol);
    validate_contains<"FC1=CC=C(C=C1)N2C(=O)C(=O)C3=CC(=CC=C23)S(=O)(=O)N4CCCC4COC5=C(F)C=C(F)C=C5">(mol);
    validate_contains<"FC1=CC=C(C=C2NC3=C(C=CC=C3)C2=N)C=C1">(mol);
    validate_contains<"FC1=CC=C(CN2C(=O)C(=O)C3=CC(=CC=C23)S(=O)(=O)N4CCCC4)C=C1">(mol);
    validate_contains<"FC1=CC=C(CN2C(=O)C(=O)C3=CC(=CC=C23)S(=O)(=O)N4CCCC4COC5=CC(F)=CC=C5)C=C1">(mol);
    validate_contains<"FC1=CC=C(CN2C(=O)C(=O)C3=CC(=CC=C23)S(=O)(=O)N4CCCC4COC5=CC=C(F)C=C5)C=C1">(mol);
    validate_contains<"FC1=CC=C(CN2C(=O)C(=O)C3=CC(=CC=C23)S(=O)(=O)N4CCCC4COC5=CC=C(F)C=C5F)C=C1">(mol);
    validate_contains<"FC1=CC=C(NC2=NC=NC3=CC=CC=C23)C=C1Cl">(mol);
    validate_contains<"FC1=CC=CC(F)=C1F">(mol);
    validate_contains<"FC1=CC=CC(NC(=O)NC2=CC=C(OC3=CC=NC4=C3NC(=O)N4)C=C2)=C1">(mol);
    validate_contains<"FC1=CN=CC=C1Cl">(mol);
    validate_contains<"FC1=NC=CN1">(mol);
    validate_contains<"FC1CCC2CCNC2C1">(mol);
    validate_contains<"FC1CCN(C1)C(=O)CNC2CCN(CC(=O)NC3=CC=C(C=C3)C#N)CC2">(mol);
    validate_contains<"FCC1=CC=CC=C1">(mol);
    validate_contains<"FCC1OCC2=C1C3=C(C=CC=C3)C=C2">(mol);
    validate_contains<"FS(=O)(=O)CC1=CC=CC=C1">(mol);
    //validate_contains<"Fc1cc(ccc1N2CCOCC2=O)N3C[C@H](CNC(=O)c4ccc(Cl)s4)OC3=O">(mol); // FIXME: stereo
    validate_contains<"Fc1ccc(Nc2ncnc3cc(OCCCN4CCOCC4)c(NC(=O)C=C)cc23)cc1Cl">(mol);
}
