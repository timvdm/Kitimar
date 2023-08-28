#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_46(Mol &mol)
{
    // SMARTS 226 - 230
    validate<"Cl(~O)(~O)(~O)(~O)">(mol);
    validate<"C(~O)(~O)(~[OH])">(mol);
    validate<"N[F,I,Br,Cl]">(mol);
    validate<"P[F,I,Br,Cl]">(mol);
    validate<"S[F,I,Br,Cl]">(mol);
}

template void RDKit_smarts_part_46<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_46<RDKit::ROMol>(RDKit::ROMol &mol);
