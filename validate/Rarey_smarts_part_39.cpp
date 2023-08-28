#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_39(Mol &mol)
{
    // SMARTS 191 - 195
    validate<"S-@P">(mol);
    validate<"S-@S">(mol);
    validate<"S=!@P">(mol);
    validate<"S=!@S">(mol);
    validate<"S=@S">(mol);
}

template void Rarey_smarts_part_39<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_39<RDKit::ROMol>(RDKit::ROMol &mol);
