#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_86(Mol &mol)
{
    // SMARTS 426 - 428
    validate<"[$([CX4,c][CH]=S),$([CX4,c]C(=S)[CX4,c])]">(mol);
    validate<"[D2R0]-[D2R0]-[D2R0]-[D2R0]-[D2R0]">(mol);
    validate<"NC(=O)N">(mol);
}

template void RDKit_smarts_part_86<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_86<RDKit::ROMol>(RDKit::ROMol &mol);
