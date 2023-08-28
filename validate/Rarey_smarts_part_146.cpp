#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_146(Mol &mol)
{
    // SMARTS 726 - 730
    validate<"[CH2X4][CX3](=[OX1])[OH0-,OH]">(mol);
    validate<"[CH2X4][OX2H]">(mol);
    validate<"[CH2X4][SX2H,SX1H0-]">(mol);
    validate<"[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1">(mol);
    validate<"[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1">(mol);
}

template void Rarey_smarts_part_146<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_146<RDKit::ROMol>(RDKit::ROMol &mol);
