#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_170(Mol &mol)
{
    // SMARTS 846 - 850
    validate<"[OD2H0][CX4][OD2H0]">(mol);
    validate<"[OD2]([#6])[#6]">(mol);
    validate<"[OH1]C1=CC(=O)NO1">(mol);
    validate<"[OH1]C1=CC(=O)ON1">(mol);
    validate<"[OH1]C1=COC=CC1=O">(mol);
}

template void Rarey_smarts_part_170<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_170<RDKit::ROMol>(RDKit::ROMol &mol);
