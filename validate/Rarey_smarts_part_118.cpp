#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_118(Mol &mol)
{
    // SMARTS 586 - 590
    validate<"[$([CX2]#C)]">(mol);
    validate<"[$([CX2](=C)=C)]">(mol);
    validate<"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]">(mol);
    validate<"[$([CX3]=[CX3])]">(mol);
    validate<"[$([CX3]=[OX1]),$([CX3+]-[OX1-])]">(mol);
}

template void Rarey_smarts_part_118<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_118<RDKit::ROMol>(RDKit::ROMol &mol);
