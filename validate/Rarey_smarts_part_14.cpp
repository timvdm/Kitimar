#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_14(Mol &mol)
{
    // SMARTS 66 - 70
    validate<"C-@P">(mol);
    validate<"C-@S">(mol);
    validate<"C-C[Br]">(mol);
    validate<"C-C[Cl]">(mol);
    validate<"C-C[F]">(mol);
}

template void Rarey_smarts_part_14<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_14<RDKit::ROMol>(RDKit::ROMol &mol);
