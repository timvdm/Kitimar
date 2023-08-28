#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_84(Mol &mol)
{
    // SMARTS 416 - 420
    validate<"[#6]S(~O)(~O)O[#6]">(mol);
    validate<"[#6][SD4](~O)(~O)[Cl,Br]">(mol);
    validate<"[ND3]([CX4])([CX4])[CX4]">(mol);
    validate<"s1c(=N)nnc1[S,N]">(mol);
    validate<"S=C([c,CX4])[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol);
}

template void RDKit_smarts_part_84<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_84<RDKit::ROMol>(RDKit::ROMol &mol);
