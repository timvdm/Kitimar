#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_37(Mol &mol)
{
    // SMARTS 181 - 185
    validate<"N(-O)(-O)">(mol);
    validate<"C1-S-C-S1">(mol);
    validate<"N=C=S">(mol);
    validate<"N=C=N">(mol);
    validate<"C1-N=N1">(mol);
}

template void RDKit_smarts_part_37<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_37<RDKit::ROMol>(RDKit::ROMol &mol);
