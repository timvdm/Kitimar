#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_32(Mol &mol)
{
    // SMARTS 156 - 160
    validate<"c-S(=O)(=O)-C=C">(mol);
    validate<"c-S(=O)(=O)-O">(mol);
    validate<"n1(=O)c([F,Br,I,Cl])cccc1">(mol);
    validate<"n1c([F,Br,I,Cl])cccc1">(mol);
    validate<"n1c([F,Br,I,Cl])ncccc1">(mol);
}

template void RDKit_smarts_part_32<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_32<RDKit::ROMol>(RDKit::ROMol &mol);
