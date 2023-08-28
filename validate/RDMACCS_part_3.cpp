#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_3(Mol &mol)
{
    // SMARTS 11 - 15
    validate<"[Cu,Zn,Ag,Cd,Au,Hg]">(mol);
    validate<"[#8]~[#7](~[#6])~[#6]">(mol);
    validate<"[#16]-[#16]">(mol);
    validate<"[#8]~[#6](~[#8])~[#8]">(mol);
    validate<"[!#6;!#1]~1~*~*~1">(mol);
}

template void RDMACCS_part_3<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_3<RDKit::ROMol>(RDKit::ROMol &mol);
