#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_24(Mol &mol)
{
    // SMARTS 116 - 120
    validate<"N(=O)-C([F,Br,I,Cl])">(mol);
    validate<"[N!H0]-[C!H0]-[N!H0]">(mol);
    validate<"N(O)=C-C-N(=O)">(mol);
    validate<"N-C(-S)=N-S(=O)(=O)">(mol);
    validate<"N-C(=N)-S">(mol);
}

template void RDKit_smarts_part_24<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_24<RDKit::ROMol>(RDKit::ROMol &mol);
