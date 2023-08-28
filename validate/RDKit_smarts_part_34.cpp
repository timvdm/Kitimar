#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_34(Mol &mol)
{
    // SMARTS 166 - 170
    validate<"C=C-N(=O)(-O)">(mol);
    validate<"C(=O)-C(=N)">(mol);
    validate<"C=C-S">(mol);
    validate<"C(-S)(-S)(-S)(-S)">(mol);
    validate<"S(-O)(-O)(-O)(-O)">(mol);
}

template void RDKit_smarts_part_34<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_34<RDKit::ROMol>(RDKit::ROMol &mol);
