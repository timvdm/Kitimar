#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_235(Mol &mol)
{
    // SMARTS 1171 - 1175
    validate<"cN=[N+]=[N-]">(mol);
    validate<"n(:c)(:c):a">(mol);
    validate<"n-!@C">(mol);
    validate<"n-!@N">(mol);
    validate<"n-!@O">(mol);
}

template void Rarey_smarts_part_235<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_235<RDKit::ROMol>(RDKit::ROMol &mol);
