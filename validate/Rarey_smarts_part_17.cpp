#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_17(Mol &mol)
{
    // SMARTS 81 - 85
    validate<"C=!@O">(mol);
    validate<"C=!@P">(mol);
    validate<"C=!@S">(mol);
    validate<"C=@A[I]">(mol);
    validate<"C=@C">(mol);
}

template void Rarey_smarts_part_17<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_17<RDKit::ROMol>(RDKit::ROMol &mol);
