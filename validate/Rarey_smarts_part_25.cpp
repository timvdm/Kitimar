#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_25(Mol &mol)
{
    // SMARTS 121 - 125
    validate<"N-@O">(mol);
    validate<"N-@P">(mol);
    validate<"N-@[Br]">(mol);
    validate<"N-@[Cl]">(mol);
    validate<"N-@[F]">(mol);
}

template void Rarey_smarts_part_25<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_25<RDKit::ROMol>(RDKit::ROMol &mol);
