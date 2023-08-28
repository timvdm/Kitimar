#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_67(Mol &mol)
{
    // SMARTS 331 - 335
    validate<"C=CC(=O)[OH]">(mol);
    validate<"N#C-C(=O)[OH]">(mol);
    validate<"N(~O)(~O)-C-C(=O)[OH]">(mol);
    validate<"[ND3]([CX4,c,H])([CX4,c,H])[CX4][$([CH]),$(C([CX4,c]))]=O">(mol);
    validate<"[OH][CX4][$([CH]),$(C([CX4,c]))]=O">(mol);
}

template void RDKit_smarts_part_67<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_67<RDKit::ROMol>(RDKit::ROMol &mol);
