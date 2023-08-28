#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_26(Mol &mol)
{
    // SMARTS 126 - 130
    validate<"N1-C-C1">(mol);
    validate<"N1-N-C(=O)-N-N1">(mol);
    validate<"N=C(-N)-C(=N)-N">(mol);
    validate<"N=C(-O)-N">(mol);
    validate<"N=C(-S)-C(=N)-N">(mol);
}

template void RDKit_smarts_part_26<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_26<RDKit::ROMol>(RDKit::ROMol &mol);
