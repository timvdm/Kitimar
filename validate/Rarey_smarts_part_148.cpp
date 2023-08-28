#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_148(Mol &mol)
{
    // SMARTS 736 - 740
    validate<"[CH2][NH2]">(mol);
    validate<"[CH3X4]">(mol);
    validate<"[CHR]=[CR][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"[CHX4]([CH3X4])[CH2X4][CH3X4]">(mol);
    validate<"[CHX4]([CH3X4])[CH3X4]">(mol);
}

template void Rarey_smarts_part_148<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_148<RDKit::ROMol>(RDKit::ROMol &mol);
