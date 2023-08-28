#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_145(Mol &mol)
{
    // SMARTS 721 - 725
    validate<"[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]">(mol);
    validate<"[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]">(mol);
    validate<"[CH2X4][CH2X4][SX2][CH3X4]">(mol);
    validate<"[CH2X4][CHX4]([CH3X4])[CH3X4]">(mol);
    validate<"[CH2X4][CX3](=[OX1])[NX3H2]">(mol);
}

template void Rarey_smarts_part_145<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_145<RDKit::ROMol>(RDKit::ROMol &mol);
