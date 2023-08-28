#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_22(Mol &mol)
{
    // SMARTS 106 - 110
    validate<"[C;!R]=[C;!R]-[S;!R]">(mol);
    validate<"C=C-S(=O)(=O)">(mol);
    validate<"C=N-C(-O)">(mol);
    validate<"C=N-C(=O)">(mol);
    validate<"C=N-N">(mol);
}

template void RDKit_smarts_part_22<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_22<RDKit::ROMol>(RDKit::ROMol &mol);
