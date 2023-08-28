#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_112(Mol &mol)
{
    // SMARTS 556 - 560
    validate<"[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]">(mol);
    validate<"[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]">(mol);
    validate<"[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]">(mol);
    validate<"[$([#16X3]=[OX1]),$([#16X3+][OX1-])]">(mol);
    validate<"[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]">(mol);
}

template void Rarey_smarts_part_112<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_112<RDKit::ROMol>(RDKit::ROMol &mol);
