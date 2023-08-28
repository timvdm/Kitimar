#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_34(Mol &mol)
{
    // SMARTS 166 - 170
    validate<"P(OCC)(OCC)(=O)C#N">(mol);
    validate<"P-!@O">(mol);
    validate<"P-!@P">(mol);
    validate<"P-!@[Br]">(mol);
    validate<"P-!@[Cl]">(mol);
}

template void Rarey_smarts_part_34<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_34<RDKit::ROMol>(RDKit::ROMol &mol);
