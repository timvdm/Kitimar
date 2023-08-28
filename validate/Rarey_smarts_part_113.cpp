#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_113(Mol &mol)
{
    // SMARTS 561 - 565
    validate<"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]">(mol);
    validate<"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]">(mol);
    validate<"[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]">(mol);
    validate<"[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]">(mol);
    validate<"[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]">(mol);
}

template void Rarey_smarts_part_113<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_113<RDKit::ROMol>(RDKit::ROMol &mol);
