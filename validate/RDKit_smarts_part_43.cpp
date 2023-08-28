#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_43(Mol &mol)
{
    // SMARTS 211 - 215
    validate<"[C;!R]=N-N=[C;!R]">(mol);
    validate<"N(C)(C)-[C;H2]-[C;H2]([F,Br,I,Cl])">(mol);
    validate<"[C;a]-[N;H1]-[N;H2]">(mol);
    validate<"[C;a]-[C;H2]([F,Br,I,Cl])">(mol);
    validate<"[S;!R]-[C;!R]-[O;!R]">(mol);
}

template void RDKit_smarts_part_43<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_43<RDKit::ROMol>(RDKit::ROMol &mol);
