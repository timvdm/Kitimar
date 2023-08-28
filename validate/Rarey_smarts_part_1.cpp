#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_1(Mol &mol)
{
    // SMARTS 1 - 5
    validate<"*!:[a]:*:*:[a]!:*">(mol);
    validate<"*!@*">(mol);
    validate<"*!@[#8]!@*">(mol);
    validate<"*(!@*)(!@*)">(mol);
    validate<"*-!:[a]:*:*:[a]-!:*">(mol);
}

template void Rarey_smarts_part_1<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_1<RDKit::ROMol>(RDKit::ROMol &mol);
