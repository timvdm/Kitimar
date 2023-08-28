#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_41(Mol &mol)
{
    // SMARTS 201 - 205
    validate<"C(=O)-N-O-C(=O)">(mol);
    validate<"C1(=O)-C=C-C(=O)-C=C1">(mol);
    validate<"B">(mol);
    validate<"ONO">(mol);
    validate<"ON(~O)~O">(mol);
}

template void RDKit_smarts_part_41<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_41<RDKit::ROMol>(RDKit::ROMol &mol);
