#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_20(Mol &mol)
{
    // SMARTS 96 - 100
    validate<"C=C-C(=O)-C=C">(mol);
    validate<"C=C-C(=S)">(mol);
    validate<"C=C-C(=S)-S">(mol);
    validate<"C=C-C([F,Br,I,Cl])=C([F,Br,I,Cl])">(mol);
    validate<"C=C-C=N">(mol);
}

template void RDKit_smarts_part_20<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_20<RDKit::ROMol>(RDKit::ROMol &mol);
