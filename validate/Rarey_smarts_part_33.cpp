#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_33(Mol &mol)
{
    // SMARTS 161 - 165
    validate<"OO">(mol);
    validate<"OS(=O)(=O)C(F)(F)F">(mol);
    validate<"P(=O)([OH])OP(=O)[OH]">(mol);
    validate<"P(=O)[O,S]">(mol);
    validate<"P(=S)(S)S">(mol);
}

template void Rarey_smarts_part_33<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_33<RDKit::ROMol>(RDKit::ROMol &mol);
