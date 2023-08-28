#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_18(Mol &mol)
{
    // SMARTS 86 - 90
    validate<"[C;H1](=O)">(mol);
    validate<"C([F,Br,I,Cl])=N">(mol);
    validate<"[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]">(mol);
    validate<"C-N=O">(mol);
    validate<"C-S(=O)(=O)-O">(mol);
}

template void RDKit_smarts_part_18<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_18<RDKit::ROMol>(RDKit::ROMol &mol);
