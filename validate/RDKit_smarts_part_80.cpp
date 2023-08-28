#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_80(Mol &mol)
{
    // SMARTS 396 - 400
    validate<"*-N(=O)(~O)">(mol);
    validate<"[N;D4]">(mol);
    validate<"OC(=O)C(=O)O">(mol);
    validate<"C=[N;R0]-[OH]">(mol);
    validate<"O~O">(mol);
}

template void RDKit_smarts_part_80<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_80<RDKit::ROMol>(RDKit::ROMol &mol);
