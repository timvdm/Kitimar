#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_164(Mol &mol)
{
    // SMARTS 816 - 820
    validate<"[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]">(mol);
    validate<"[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]">(mol);
    validate<"[NX3;H2,H1;!$(NC=O)]">(mol);
    validate<"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]">(mol);
    validate<"[NX3H2,NH3X4+][CX4H]([*])[CX3](=[OX1])[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-]">(mol);
}

template void Rarey_smarts_part_164<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_164<RDKit::ROMol>(RDKit::ROMol &mol);
