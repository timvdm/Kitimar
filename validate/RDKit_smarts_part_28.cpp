#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_28(Mol &mol)
{
    // SMARTS 136 - 140
    validate<"N=C=S">(mol);
    validate<"O-C(=O)-O-N">(mol);
    validate<"O-C([F,Br,I,Cl])=S">(mol);
    validate<"O-C=C">(mol);
    validate<"O-N=C-C=N-O">(mol);
}

template void RDKit_smarts_part_28<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_28<RDKit::ROMol>(RDKit::ROMol &mol);
