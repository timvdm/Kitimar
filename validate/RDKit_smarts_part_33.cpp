#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_33(Mol &mol)
{
    // SMARTS 161 - 165
    validate<"C1-C-O-C-O-C1">(mol);
    validate<"C=S">(mol);
    validate<"c=S">(mol);
    validate<"S=C-S">(mol);
    validate<"C=C-N(=O)(=O)">(mol);
}

template void RDKit_smarts_part_33<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_33<RDKit::ROMol>(RDKit::ROMol &mol);
