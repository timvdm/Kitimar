#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_28(Mol &mol)
{
    // SMARTS 136 - 140
    validate<"N=NC(S)N">(mol);
    validate<"N=[O,C,N,S]">(mol);
    validate<"NP(=O)(N)N">(mol);
    validate<"N[CH2]C#N">(mol);
    validate<"N[CX4H2][CX3](=[OX1])[O,N]">(mol);
}

template void Rarey_smarts_part_28<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_28<RDKit::ROMol>(RDKit::ROMol &mol);
