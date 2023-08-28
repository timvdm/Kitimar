#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_60(Mol &mol)
{
    // SMARTS 296 - 300
    validate<"[C;R0;X4]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH][CH2]">(mol);
    validate<"[CH]=O">(mol);
    validate<"[C;R0;X4]!@[CH2]!@[CH2]!@[CH2]!@[CH2]!@[CH2]!@[C;R0;X4]">(mol);
    validate<"[S;D2;R0]-[S;D2]">(mol);
    validate<"[SH]">(mol);
}

template void RDKit_smarts_part_60<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_60<RDKit::ROMol>(RDKit::ROMol &mol);
