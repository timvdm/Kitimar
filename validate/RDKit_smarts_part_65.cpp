#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_65(Mol &mol)
{
    // SMARTS 321 - 325
    validate<"[CX4;H][OH]">(mol);
    validate<"[CX4;H0][OH]">(mol);
    validate<"[#6][CH]=O">(mol);
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]">(mol);
    validate<"[CX4][Cl,Br,I]">(mol);
}

template void RDKit_smarts_part_65<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_65<RDKit::ROMol>(RDKit::ROMol &mol);
