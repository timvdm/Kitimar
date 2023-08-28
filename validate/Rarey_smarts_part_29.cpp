#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_29(Mol &mol)
{
    // SMARTS 141 - 145
    validate<"N~C(~O)~N">(mol);
    validate<"N~C~C(~O)~O">(mol);
    validate<"O!@[C$(C@*)]">(mol);
    validate<"O-!@O">(mol);
    validate<"O-!@S">(mol);
}

template void Rarey_smarts_part_29<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_29<RDKit::ROMol>(RDKit::ROMol &mol);
