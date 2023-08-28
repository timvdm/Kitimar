#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_15(Mol &mol)
{
    // SMARTS 71 - 75
    validate<"C(=O)-C(-O)-N-C=O">(mol);
    validate<"C(=O)-C(=C)-C(=O)">(mol);
    validate<"C(=O)-C(=N)">(mol);
    validate<"C(=O)-C([F,Br,I,Cl])-C(=O)">(mol);
    validate<"C(=O)-C([F,Br,I,Cl])=C([F,Br,I,Cl])-C(=O)">(mol);
}

template void RDKit_smarts_part_15<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_15<RDKit::ROMol>(RDKit::ROMol &mol);
