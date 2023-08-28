#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_29(Mol &mol)
{
    // SMARTS 141 - 145
    validate<"[!#1]!:*:*!:[!#1]">(mol);
    validate<"*~1~*~*~*~*~*~1">(mol);
    validate<"[#8]">(mol);
    validate<"[$([!#1]~[CH2]~[CH2]~[!#1]),$([R]1@[CH2;R]@[CH2;R]1)]">(mol);
    validate<"[!#1]~[!#6;!#1](~[!#1])~[!#1]">(mol);
}

template void RDMACCS_part_29<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_29<RDKit::ROMol>(RDKit::ROMol &mol);
