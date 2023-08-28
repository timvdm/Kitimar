#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_23(Mol &mol)
{
    // SMARTS 111 - 115
    validate<"[#8]!:*:*">(mol);
    validate<"[CH3]~[CH2]~[!#1]">(mol);
    validate<"[CH3]~*~[CH2]~[!#1]">(mol);
    validate<"[$([CH3]~*~*~[CH2]~[!#1]),$([CH3]~*1~*~[CH2]1)]">(mol);
    validate<"[#7]~*~[#8]">(mol);
}

template void RDMACCS_part_23<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_23<RDKit::ROMol>(RDKit::ROMol &mol);
