#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_21(Mol &mol)
{
    // SMARTS 101 - 105
    validate<"C=C-N([C,c])([C,c])">(mol);
    validate<"C=C-N(=O)(=O)">(mol);
    validate<"C=C-N(O)(=O)">(mol);
    validate<"C=C-O">(mol);
    validate<"C=C-O-C(=O)">(mol);
}

template void RDKit_smarts_part_21<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_21<RDKit::ROMol>(RDKit::ROMol &mol);
