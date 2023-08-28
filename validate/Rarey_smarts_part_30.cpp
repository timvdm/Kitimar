#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_30(Mol &mol)
{
    // SMARTS 146 - 150
    validate<"O-!@[Br]">(mol);
    validate<"O-!@[Cl]">(mol);
    validate<"O-!@[F]">(mol);
    validate<"O-!@[I]">(mol);
    validate<"O-@O">(mol);
}

template void Rarey_smarts_part_30<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_30<RDKit::ROMol>(RDKit::ROMol &mol);
