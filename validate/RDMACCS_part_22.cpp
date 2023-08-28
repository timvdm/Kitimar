#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_22(Mol &mol)
{
    // SMARTS 106 - 110
    validate<"[CH3]~*~*~*~[CH2]~[!#1]">(mol);
    validate<"[!#1]~[CH2]~[#8]">(mol);
    validate<"[#7]~[#6]~[#8]">(mol);
    validate<"[#7]~*~[CH2]~[!#1]">(mol);
    validate<"[!#1]~*(~[!#1])(~[!#1])~[!#1]">(mol);
}

template void RDMACCS_part_22<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_22<RDKit::ROMol>(RDKit::ROMol &mol);
