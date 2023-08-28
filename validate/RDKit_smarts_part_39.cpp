#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_39(Mol &mol)
{
    // SMARTS 191 - 195
    validate<"C=C-N">(mol);
    validate<"C=N-S(=O)(=O)">(mol);
    validate<"N-N=O">(mol);
    validate<"N-C(=O)-S">(mol);
    validate<"S-C#N">(mol);
}

template void RDKit_smarts_part_39<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_39<RDKit::ROMol>(RDKit::ROMol &mol);
