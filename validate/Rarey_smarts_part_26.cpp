#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_26(Mol &mol)
{
    // SMARTS 126 - 130
    validate<"N-@[I]">(mol);
    validate<"N-[C,c,N]=[C,c,N,n,O,S]">(mol);
    validate<"N1CCC1=O">(mol);
    validate<"N1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[O,N]">(mol);
    validate<"N=!@N">(mol);
}

template void Rarey_smarts_part_26<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_26<RDKit::ROMol>(RDKit::ROMol &mol);
