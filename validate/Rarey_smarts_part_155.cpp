#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_155(Mol &mol)
{
    // SMARTS 771 - 775
    validate<"[HC]=O">(mol);
    validate<"[H]">(mol);
    validate<"[H][H]">(mol);
    validate<"[N&D2](=O)">(mol);
    validate<"[N&X4&+,N&X3&+0]">(mol);
}

template void Rarey_smarts_part_155<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_155<RDKit::ROMol>(RDKit::ROMol &mol);
