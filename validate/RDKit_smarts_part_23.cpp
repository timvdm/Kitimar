#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_23(Mol &mol)
{
    // SMARTS 111 - 115
    validate<"C=N-N=C">(mol);
    validate<"C=N-O">(mol);
    validate<"C=N-O-C(=O)">(mol);
    validate<"C=N=O">(mol);
    validate<"N#C-S">(mol);
}

template void RDKit_smarts_part_23<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_23<RDKit::ROMol>(RDKit::ROMol &mol);
