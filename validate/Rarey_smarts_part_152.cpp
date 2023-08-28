#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_152(Mol &mol)
{
    // SMARTS 756 - 760
    validate<"[CX3](=[OX1])[NX3H][CX3](=[OX1])">(mol);
    validate<"[CX3](=[OX1])[OX2][CX3](=[OX1])">(mol);
    validate<"[CX3]=[CD2H1]-[CD2H1]=[O,N,S]">(mol);
    validate<"[CX3]=[OX1]">(mol);
    validate<"[CX4]">(mol);
}

template void Rarey_smarts_part_152<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_152<RDKit::ROMol>(RDKit::ROMol &mol);
