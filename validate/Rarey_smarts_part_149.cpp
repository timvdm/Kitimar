#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_149(Mol &mol)
{
    // SMARTS 741 - 745
    validate<"[CHX4]([CH3X4])[OX2H]">(mol);
    validate<"[CH]#C">(mol);
    validate<"[CX1-]#[NX2+]">(mol);
    validate<"[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]">(mol);
    validate<"[CX3H1](=O)[#6]">(mol);
}

template void Rarey_smarts_part_149<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_149<RDKit::ROMol>(RDKit::ROMol &mol);
