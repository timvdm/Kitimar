#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_45(Mol &mol)
{
    // SMARTS 221 - 225
    validate<"[N+]#N">(mol);
    validate<"ON#C">(mol);
    validate<"OOO">(mol);
    validate<"OO">(mol);
    validate<"Cl(~O)(~O)(~O)">(mol);
}

template void RDKit_smarts_part_45<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_45<RDKit::ROMol>(RDKit::ROMol &mol);
