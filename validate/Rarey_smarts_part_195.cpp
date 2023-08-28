#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_195(Mol &mol)
{
    // SMARTS 971 - 975
    validate<"c-!@O">(mol);
    validate<"c-!@S">(mol);
    validate<"c-@C">(mol);
    validate<"c-@N">(mol);
    validate<"c-@O">(mol);
}

template void Rarey_smarts_part_195<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_195<RDKit::ROMol>(RDKit::ROMol &mol);
