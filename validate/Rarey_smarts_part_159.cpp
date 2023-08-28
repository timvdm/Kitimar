#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_159(Mol &mol)
{
    // SMARTS 791 - 795
    validate<"[ND2H0]#[CD1H0]">(mol);
    validate<"[ND2H0]=[CD2H0]=[ND2H0]">(mol);
    validate<"[ND2H0]=[CD2H0]=[OD1H0,SD1H0]">(mol);
    validate<"[ND2H0]=[ND1H0]">(mol);
    validate<"[ND2H0]=[ND2H0]">(mol);
}

template void Rarey_smarts_part_159<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_159<RDKit::ROMol>(RDKit::ROMol &mol);
