#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_114(Mol &mol)
{
    // SMARTS 566 - 570
    validate<"[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]">(mol);
    validate<"[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]">(mol);
    validate<"[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]">(mol);
    validate<"[$([#1X1][$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)])]">(mol);
    validate<"[$([#1X1][$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)])]">(mol);
}

template void Rarey_smarts_part_114<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_114<RDKit::ROMol>(RDKit::ROMol &mol);
