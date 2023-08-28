#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_25(Mol &mol)
{
    // SMARTS 121 - 125
    validate<"N-C(=O)-O-c1ccccc1">(mol);
    validate<"N-C(=S)">(mol);
    validate<"N-C(=S)-N">(mol);
    validate<"N-C=C">(mol);
    validate<"N-O-C(=O)">(mol);
}

template void RDKit_smarts_part_25<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_25<RDKit::ROMol>(RDKit::ROMol &mol);
