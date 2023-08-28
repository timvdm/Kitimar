#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_35(Mol &mol)
{
    // SMARTS 171 - 175
    validate<"S(=O)(=O)-C-N(=O)(=O)">(mol);
    validate<"S(=O)(=O)-C-N(=O)(-O)">(mol);
    validate<"N#C-S">(mol);
    validate<"C=N-O">(mol);
    validate<"C=N-S=N">(mol);
}

template void RDKit_smarts_part_35<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_35<RDKit::ROMol>(RDKit::ROMol &mol);
