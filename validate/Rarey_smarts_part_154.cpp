#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_154(Mol &mol)
{
    // SMARTS 766 - 770
    validate<"[Cl]C([C&R0])=N">(mol);
    validate<"[F,Cl,Br,I]">(mol);
    validate<"[H+]">(mol);
    validate<"[H,#1]">(mol);
    validate<"[H1,H0-]">(mol);
}

template void Rarey_smarts_part_154<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_154<RDKit::ROMol>(RDKit::ROMol &mol);
