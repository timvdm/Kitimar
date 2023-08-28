#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_147(Mol &mol)
{
    // SMARTS 731 - 735
    validate<"[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12">(mol);
    validate<"[CH2]=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"[CH2]=O">(mol);
    validate<"[CH2]N([CH3])[CH3]">(mol);
    validate<"[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]">(mol);
}

template void Rarey_smarts_part_147<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_147<RDKit::ROMol>(RDKit::ROMol &mol);
