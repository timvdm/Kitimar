#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_17(Mol &mol)
{
    // SMARTS 81 - 85
    validate<"C(=O)-N-N-C(=O)">(mol);
    validate<"C(=O)-O-N-C(=O)">(mol);
    validate<"C(=O)-S">(mol);
    validate<"C(=O)-S-C(=S)">(mol);
    validate<"C(=S)-S">(mol);
}

template void RDKit_smarts_part_17<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_17<RDKit::ROMol>(RDKit::ROMol &mol);
