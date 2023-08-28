#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_50(Mol &mol)
{
    // SMARTS 246 - 250
    validate<"[C;H2]=C-C(=O)">(mol);
    validate<"[C;D1&H3,D2&H2,D3&H1,D4]-[C;H1]=C-C(=O)-C">(mol);
    validate<"[C;D1&H3,D2&H2,D3&H1,D4]-[C;H1]=C-C(=O)-O-C-C">(mol);
    validate<"C#C-C#N">(mol);
    validate<"C#C-C(=O)">(mol);
}

template void RDKit_smarts_part_50<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_50<RDKit::ROMol>(RDKit::ROMol &mol);
