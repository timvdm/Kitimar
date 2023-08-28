#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_4(Mol &mol)
{
    // SMARTS 16 - 20
    validate<"*1*******1">(mol);
    validate<"*1******1">(mol);
    validate<"*1*****1">(mol);
    validate<"*1****1">(mol);
    validate<"*1***1">(mol);
}

template void Rarey_smarts_part_4<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_4<RDKit::ROMol>(RDKit::ROMol &mol);
