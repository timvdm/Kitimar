#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_35(Mol &mol)
{
    // SMARTS 171 - 175
    validate<"P-!@[F]">(mol);
    validate<"P-!@[I]">(mol);
    validate<"P-@O">(mol);
    validate<"P-@P">(mol);
    validate<"P=!@P">(mol);
}

template void Rarey_smarts_part_35<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_35<RDKit::ROMol>(RDKit::ROMol &mol);
