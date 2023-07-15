#include "Validate.hpp"

void BindingDB_substructure_part_35(OpenBabel::OBMol &mol)
{
    // SMARTS 1701 - 1750
    validate_contains<"COC1=CC(Br)=CC(C(=O)NC2=CC=CC=C2)=C1OC">(mol);
    //validate_contains<"COC1=CC(C=C2/CCC(=C/C3=CC(OC)=C(O)C=C3)C2=O)=CC=C1O">(mol); // FIXME: stereo
    validate_contains<"COC1=CC(C=CC)=CC=C1O">(mol);
    validate_contains<"COC1=CC(C=O)=CC=C1O">(mol);
    validate_contains<"COC1=CC(CC2=CN=C(N)N=C2N)=CC(OC)=C1OC">(mol);
    validate_contains<"COC1=CC(CN)=C(OC)C=C1">(mol);
    validate_contains<"COC1=CC(Cl)=CC(C(=O)NC2=NC=C(Cl)C=C2)=C1NC(=O)C3=C(Cl)C(CN(C)C4=NCCO4)=CS3">(mol);
    validate_contains<"COC1=CC(NC2=C(C=NC3=CC(OCCCN4CCN(C)CC4)=CC=C23)C#N)=C(Cl)C=C1Cl">(mol);
    validate_contains<"COC1=CC(O)=CC(OC)=C1OC">(mol);
    validate_contains<"COC1=CC2=C(C(=O)C(OC)=C(O2)C3=CC=C(O)C=C3)C(O)=C1">(mol);
    validate_contains<"COC1=CC2=C(C=C1)N=C(S2)C3=CC(Br)=C(N)C=C3">(mol);
    validate_contains<"COC1=CC2=C(C=C1)N=CC=C2O">(mol);
    validate_contains<"COC1=CC2=C(C=C1OC)C(=O)C(C)C(C)C2C3=CC(OC)=C(OC)C=C3">(mol);
    validate_contains<"COC1=CC2=C(C=C1OC)C3=NN(CC3C2C4=CC=C(Cl)C=C4)C(=O)C5=CC=NC=C5">(mol);
    validate_contains<"COC1=CC2=C(N)C3=C(C=C(Cl)C=C3)N=C2C=C1">(mol);
    validate_contains<"COC1=CC2=C(NC=C2)C=C1">(mol);
    validate_contains<"COC1=CC=C(C=C1)C(=O)C2=CC(OC)=C(OC)C(OC)=C2">(mol);
    validate_contains<"COC1=CC=C(C=C1)C(=O)CC2=CC=C(O)C=C2">(mol);
    validate_contains<"COC1=CC=C(C=C1)C(=O)NC2=CC=C(O)C=C2NC(=O)C3=CC=C(C=C3)N(C)C">(mol);
    validate_contains<"COC1=CC=C(C=C1)C(CN(C)C)C2(O)CCCCC2">(mol);
    validate_contains<"COC1=CC=C(C=C1)C(NC(=O)CC2=CC(OC)=C(C=C2Br)C(O)=O)C3=C(C=CC=C3C)N4CCCCC4">(mol);
    //validate_contains<"COC1=CC=C(C=C1)C2=C(CCNS(=O)(=O)N[C@@H]3CCN(CC4=CC=CC=C4)C3)N=C(N2)C5=CC=CC=C5">(mol); // FIXME: stereo
    validate_contains<"COC1=CC=C(C=C1)C2=CC=C(C=C2)C3CCC(=N3)C4=C(F)C=CC=C4F">(mol);
    validate_contains<"COC1=CC=C(C=C1)C2=NC(=CS2)C3=CC(C(=O)NN4CCCCC4)=C(C)N3CC5CCCCC5">(mol);
    //validate_contains<"COC1=CC=C(C=C1)C=C2/CCCC(=C/C3=CC=C(OC)C=C3)C2=O">(mol); // FIXME: stereo
    validate_contains<"COC1=CC=C(C=C1)N2C=CC=N2">(mol);
    validate_contains<"COC1=CC=C(C=C1)N2N=NN=C2SCC(=O)NN=CC3=CNN=C3C4=CC=C(F)C=C4">(mol);
    validate_contains<"COC1=CC=C(C=C1)S(=O)(=O)C(CC2=CC=CC=N2)C(C(C)C)C(=O)NO">(mol);
    validate_contains<"COC1=CC=C(C=C1OC)C(=C)CC(C)C">(mol);
    validate_contains<"COC1=CC=C(C=C1OC)C2OC3=CC(O)=CC(O)=C3C(=O)C2O">(mol);
    //validate_contains<"COC1=CC=C(C=C2/CCCC2=O)C=C1">(mol); // FIXME: stereo
    validate_contains<"COC1=CC=C(CCOC2=CC=C(NC(=O)NC3=CC(=NN3C3=CC=C(C)C=C3)C(C)(C)C)C3=C2C=CC=C3)C=C1OC">(mol);
    validate_contains<"COC1=CC=C(CN)C=C1">(mol);
    validate_contains<"COC1=CC=C(CNC(=O)NC2=NC=C(S2)N(=O)=O)C=C1">(mol);
    validate_contains<"COC1=CC=C(CNC(N)=O)C=C1">(mol);
    validate_contains<"COC1=CC=C2C3N(CC(C4=CC=C(Cl)C=C4)C5=C3NC6=CC=CC=C56)C(=O)C2=C1OC">(mol);
    validate_contains<"COC1=CC=C2C=C(C(=O)NCC3=CC=CC=C3)C(=O)OC2=C1">(mol);
    validate_contains<"COC1=CC=C2C=C3C4=CC5=C(OCO5)C=C4CC[N+]3=CC2=C1OC">(mol);
    validate_contains<"COC1=CC=C2N=CC(C#N)=C(NC3=CC=CC=C3)C2=C1">(mol);
    validate_contains<"COC1=CC=C2N=CN=C(NC3=C(F)C=C(Br)C=C3)C2=C1">(mol);
    validate_contains<"COC1=CC=C2N=CN=C(NC3=CC=CC=C3F)C2=C1">(mol);
    validate_contains<"COC1=CC=CC(=C1)COC2=NC(=NC3=C2NC=N3)N">(mol);
    validate_contains<"COC1=CC=CC(NS(=O)(=O)C2=CC=C(N)C=C2)=C1OC">(mol);
    validate_contains<"COC1=CC=CC(OC)=C1OC">(mol);
    //validate_contains<"COC1=CC=CC2=C1C(=O)C3=C(C(O)=C4C[C@](O)(C[C@H](O[C@H]5C[C@H](N)[C@H](O)[C@H](C)O5)C4=C3O)C(=O)CO)C2=O">(mol); // FIXME: stereo
    validate_contains<"COC1=CC=CC2=C1N=CC(C)=C2">(mol);
    validate_contains<"COC1=CC=CC2=C1OC(=O)C(=C2)C(=O)NC3=CC=CC=C3">(mol);
    validate_contains<"COC1=CC=CC=C1">(mol);
    validate_contains<"COC1=CC=CC=C1C">(mol);
    validate_contains<"COC1=CC=CC=C1OC2=C(OCCO)N=C(N=C2NS(=O)(=O)C3=CC=C(C=C3)C(C)(C)C)C4=NC=CC=N4">(mol);
}
