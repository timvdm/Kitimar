#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_119(Mol &mol)
{
    // SMARTS 591 - 595
    validate<"[$([CX3]=[OX1])]">(mol);
    validate<"[$([NH2]!:c),$([NH1]([CX4])!:c),$([NH0]([CX4])([CX4])!:c)]">(mol);
    validate<"[$([NH2]-C),$([OH]-*)]">(mol);
    validate<"[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]">(mol);
    validate<"[$([NX1]#*)]">(mol);
}

template void Rarey_smarts_part_119<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_119<RDKit::ROMol>(RDKit::ROMol &mol);
