#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_2(Mol &mol)
{
    // SMARTS 6 - 10
    validate<"C(=O)OC(=O)">(mol);
    validate<"OO">(mol);
    validate<"C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)">(mol);
    validate<"C(=O)Oc1ccc(N(=O)=O)cc1">(mol);
    validate<"C(=O)Onnn">(mol);
}

template void RDKit_smarts_part_2<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_2<RDKit::ROMol>(RDKit::ROMol &mol);
