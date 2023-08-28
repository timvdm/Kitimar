#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_64(Mol &mol)
{
    // SMARTS 316 - 320
    validate<"C(=O)Cl">(mol);
    validate<"[$(C#N),$([C,N,S]=O)][CH2,CH][$([C,N,S]=O),$(C#N)]">(mol);
    validate<"C(=O)[Cl,Br]">(mol);
    validate<"[CX4][OH]">(mol);
    validate<"[CX4;H2][OH]">(mol);
}

template void RDKit_smarts_part_64<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_64<RDKit::ROMol>(RDKit::ROMol &mol);
