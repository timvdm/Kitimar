#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_6(Mol &mol)
{
    // SMARTS 26 - 30
    validate<"N#CC[OH]">(mol);
    validate<"N#CC(=O)">(mol);
    validate<"S(=O)(=O)C#N">(mol);
    validate<"P(OCC)(OCC)(=O)C#N">(mol);
    validate<"[N;R0]=[N;R0]C#N">(mol);
}

template void RDKit_smarts_part_6<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_6<RDKit::ROMol>(RDKit::ROMol &mol);
