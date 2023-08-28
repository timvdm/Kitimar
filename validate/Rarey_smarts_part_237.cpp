#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_237(Mol &mol)
{
    // SMARTS 1181 - 1185
    validate<"n-@S">(mol);
    validate<"n-[Br]">(mol);
    validate<"n-[Cl]">(mol);
    validate<"n-[F]">(mol);
    validate<"n-[I]">(mol);
}

template void Rarey_smarts_part_237<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_237<RDKit::ROMol>(RDKit::ROMol &mol);
