#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_55(Mol &mol)
{
    // SMARTS 271 - 275
    validate<"S=C=N">(mol);
    validate<"O~O">(mol);
    validate<"[Si]~N">(mol);
    validate<"P[S,N]">(mol);
    validate<"[N;R0]=[N;R0]=[C;R0]">(mol);
}

template void RDKit_smarts_part_55<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_55<RDKit::ROMol>(RDKit::ROMol &mol);
