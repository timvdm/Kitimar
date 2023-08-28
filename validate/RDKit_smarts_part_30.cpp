#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_30(Mol &mol)
{
    // SMARTS 146 - 150
    validate<"S-C(=N)-N">(mol);
    validate<"S-C(=N)-S">(mol);
    validate<"S-C(=S)-N">(mol);
    validate<"S-N-C(=O)">(mol);
    validate<"[S;H]">(mol);
}

template void RDKit_smarts_part_30<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_30<RDKit::ROMol>(RDKit::ROMol &mol);
