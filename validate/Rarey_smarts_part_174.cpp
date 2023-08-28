#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_174(Mol &mol)
{
    // SMARTS 866 - 870
    validate<"[OH1][P,C,S](=O)">(mol);
    validate<"[OH1]c1c[c,n]ccc1">(mol);
    validate<"[OH1]c1oncc1">(mol);
    validate<"[OH]c1cccc2ccccc12">(mol);
    validate<"[OX1]=CN">(mol);
}

template void Rarey_smarts_part_174<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_174<RDKit::ROMol>(RDKit::ROMol &mol);
