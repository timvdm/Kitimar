#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_236(Mol &mol)
{
    // SMARTS 1176 - 1180
    validate<"n-!@S">(mol);
    validate<"n-@C">(mol);
    validate<"n-@N">(mol);
    validate<"n-@O">(mol);
    validate<"n-@P">(mol);
}

template void Rarey_smarts_part_236<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_236<RDKit::ROMol>(RDKit::ROMol &mol);
