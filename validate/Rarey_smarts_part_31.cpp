#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_31(Mol &mol)
{
    // SMARTS 151 - 155
    validate<"O-@S">(mol);
    validate<"O1CCCCC1C2CCCO2">(mol);
    validate<"O1CCCCC1OC2CCC3CCCCC3C2">(mol);
    validate<"O=!@P">(mol);
    validate<"O=!@S">(mol);
}

template void Rarey_smarts_part_31<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_31<RDKit::ROMol>(RDKit::ROMol &mol);
