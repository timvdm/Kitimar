#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_7(Mol &mol)
{
    // SMARTS 31 - 35
    validate<"[N;R0]=[N;R0]CC=O">(mol);
    validate<"[CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0]">(mol);
    validate<"[O;R1][C;R1][C;R1][O;R1][C;R1][C;R1][O;R1]">(mol);
    validate<"SS">(mol);
    validate<"[SH]">(mol);
}

template void RDKit_smarts_part_7<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_7<RDKit::ROMol>(RDKit::ROMol &mol);
