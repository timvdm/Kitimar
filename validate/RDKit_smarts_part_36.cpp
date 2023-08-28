#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_36(Mol &mol)
{
    // SMARTS 176 - 180
    validate<"S-C(=O)-S">(mol);
    validate<"S-C(=N)-S">(mol);
    validate<"S-C(=N)(=N)">(mol);
    validate<"N(=O)(=O)">(mol);
    validate<"N(=O)(-O)">(mol);
}

template void RDKit_smarts_part_36<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_36<RDKit::ROMol>(RDKit::ROMol &mol);
