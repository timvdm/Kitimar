#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_27(Mol &mol)
{
    // SMARTS 131 - 135
    validate<"N=C([F,Br,I,Cl])">(mol);
    validate<"N=C-C(=O)">(mol);
    validate<"N=C-S">(mol);
    validate<"N=C=N">(mol);
    validate<"N=C=O">(mol);
}

template void RDKit_smarts_part_27<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_27<RDKit::ROMol>(RDKit::ROMol &mol);
