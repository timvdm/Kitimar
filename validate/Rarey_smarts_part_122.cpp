#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_122(Mol &mol)
{
    // SMARTS 606 - 610
    validate<"[$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)]">(mol);
    validate<"[$([NX4+]),$([NX4]=*)]">(mol);
    validate<"[$([OH]-*=[!#6])]">(mol);
    validate<"[$([OX1]=[CX3])]">(mol);
    validate<"[$([OX1]=[NX3](=[OX1])[OX1-]),$([OX1]=[NX3+]([OX1-])[OX1-])]">(mol);
}

template void Rarey_smarts_part_122<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_122<RDKit::ROMol>(RDKit::ROMol &mol);
