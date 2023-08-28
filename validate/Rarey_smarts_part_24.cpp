#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_24(Mol &mol)
{
    // SMARTS 116 - 120
    validate<"N-!@[Br]">(mol);
    validate<"N-!@[Cl]">(mol);
    validate<"N-!@[F]">(mol);
    validate<"N-!@[I]">(mol);
    validate<"N-@N">(mol);
}

template void Rarey_smarts_part_24<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_24<RDKit::ROMol>(RDKit::ROMol &mol);
