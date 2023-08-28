#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_151(Mol &mol)
{
    // SMARTS 751 - 755
    validate<"[CX3](=[OX1])C">(mol);
    validate<"[CX3](=[OX1])O">(mol);
    validate<"[CX3](=[OX1])[F,Cl,Br,I]">(mol);
    validate<"[CX3](=[OX1])[NX3H0]([#6])[CX3](=[OX1])">(mol);
    validate<"[CX3](=[OX1])[NX3H0]([NX3H0]([CX3](=[OX1]))[CX3](=[OX1]))[CX3](=[OX1])">(mol);
}

template void Rarey_smarts_part_151<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_151<RDKit::ROMol>(RDKit::ROMol &mol);
