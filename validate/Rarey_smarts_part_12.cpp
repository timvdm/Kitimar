#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_12(Mol &mol)
{
    // SMARTS 56 - 60
    validate<"C-!@[Cl]">(mol);
    validate<"C-!@[F]">(mol);
    validate<"C-!@[I]">(mol);
    validate<"C-@A[Br]">(mol);
    validate<"C-@A[Cl]">(mol);
}

template void Rarey_smarts_part_12<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_12<RDKit::ROMol>(RDKit::ROMol &mol);
