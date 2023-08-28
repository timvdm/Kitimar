#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_3(Mol &mol)
{
    // SMARTS 11 - 15
    validate<"N=C=[S,O]">(mol);
    validate<"OS(=O)(=O)C(F)(F)F">(mol);
    validate<"P(=S)(S)S">(mol);
    validate<"NP(=O)(N)N">(mol);
    validate<"cN=[N+]=[N-]">(mol);
}

template void RDKit_smarts_part_3<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_3<RDKit::ROMol>(RDKit::ROMol &mol);
