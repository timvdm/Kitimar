#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_9(Mol &mol)
{
    // SMARTS 41 - 45
    validate<"C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)">(mol);
    validate<"C(=O)Oc1ccc(N(=O)=O)cc1">(mol);
    validate<"C(=O)Onnn">(mol);
    validate<"C(=O)[Cl,Br,I]">(mol);
    validate<"C(=O)[OH1]">(mol);
}

template void Rarey_smarts_part_9<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_9<RDKit::ROMol>(RDKit::ROMol &mol);
