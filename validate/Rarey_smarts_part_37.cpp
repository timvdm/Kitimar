#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_37(Mol &mol)
{
    // SMARTS 181 - 185
    validate<"P[OH1]">(mol);
    validate<"S(=O)(=O)C#N">(mol);
    validate<"S(=O)(=O)[F,Cl,Br,I]">(mol);
    validate<"S([#6])[CX3](=O)[#6]">(mol);
    validate<"S-!@P">(mol);
}

template void Rarey_smarts_part_37<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_37<RDKit::ROMol>(RDKit::ROMol &mol);
