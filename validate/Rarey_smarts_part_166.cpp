#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_166(Mol &mol)
{
    // SMARTS 826 - 830
    validate<"[NX3][CX3](=[OX1])[#6]">(mol);
    validate<"[NX3][CX3](=[OX1])[OX2H0]">(mol);
    validate<"[NX3][CX3]=[CX3]">(mol);
    validate<"[NX3][CX3]=[NX3+]">(mol);
    validate<"[NX3][CX3]=[SX1]">(mol);
}

template void Rarey_smarts_part_166<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_166<RDKit::ROMol>(RDKit::ROMol &mol);
