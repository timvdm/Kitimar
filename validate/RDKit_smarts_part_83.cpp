#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_83(Mol &mol)
{
    // SMARTS 411 - 415
    validate<"C14~*~*~*~*~C~1~*~*~C2~C3~*~*~*~C~3~*~*~C~2~4">(mol);
    validate<"[#6][SD3](~O)[#6]">(mol);
    validate<"[#6][SD4](~O)(~O)[#6]">(mol);
    validate<"[#6][SD4](~O)(~O)N">(mol);
    validate<"[#6]S(~O)(~O)[OH]">(mol);
}

template void RDKit_smarts_part_83<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_83<RDKit::ROMol>(RDKit::ROMol &mol);
