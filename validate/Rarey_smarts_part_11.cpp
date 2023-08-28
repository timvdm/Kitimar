#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_11(Mol &mol)
{
    // SMARTS 51 - 55
    validate<"C-!@N">(mol);
    validate<"C-!@O">(mol);
    validate<"C-!@P">(mol);
    validate<"C-!@S">(mol);
    validate<"C-!@[Br]">(mol);
}

template void Rarey_smarts_part_11<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_11<RDKit::ROMol>(RDKit::ROMol &mol);
