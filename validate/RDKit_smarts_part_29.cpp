#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_29(Mol &mol)
{
    // SMARTS 141 - 145
    validate<"P">(mol);
    validate<"S(=O)(=O)-C([F,Br,I,Cl])">(mol);
    validate<"S(=O)(=O)O">(mol);
    validate<"S-C#N">(mol);
    validate<"S-C(=C)-S">(mol);
}

template void RDKit_smarts_part_29<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_29<RDKit::ROMol>(RDKit::ROMol &mol);
