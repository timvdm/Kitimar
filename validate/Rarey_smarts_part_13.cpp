#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_13(Mol &mol)
{
    // SMARTS 61 - 65
    validate<"C-@A[F]">(mol);
    validate<"C-@A[I]">(mol);
    validate<"C-@C">(mol);
    validate<"C-@N">(mol);
    validate<"C-@O">(mol);
}

template void Rarey_smarts_part_13<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_13<RDKit::ROMol>(RDKit::ROMol &mol);
