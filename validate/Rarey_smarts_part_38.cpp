#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_38(Mol &mol)
{
    // SMARTS 186 - 190
    validate<"S-!@S">(mol);
    validate<"S-!@[Br]">(mol);
    validate<"S-!@[Cl]">(mol);
    validate<"S-!@[F]">(mol);
    validate<"S-!@[I]">(mol);
}

template void Rarey_smarts_part_38<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_38<RDKit::ROMol>(RDKit::ROMol &mol);
