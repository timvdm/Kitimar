#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_36(Mol &mol)
{
    // SMARTS 176 - 180
    validate<"P=@P">(mol);
    validate<"P=P@[Br]">(mol);
    validate<"P=P@[Cl]">(mol);
    validate<"P=P@[F]">(mol);
    validate<"P=P@[I]">(mol);
}

template void Rarey_smarts_part_36<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_36<RDKit::ROMol>(RDKit::ROMol &mol);
