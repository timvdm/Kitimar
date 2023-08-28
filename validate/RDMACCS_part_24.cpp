#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_24(Mol &mol)
{
    // SMARTS 116 - 120
    validate<"[$([!#1]~[CH2]~[CH2]~[!#1]),$(*1~[CH2]~[CH2]1)]">(mol);
    validate<"[#7]=*">(mol);
    validate<"[!#6;R]">(mol);
    validate<"[#7;R]">(mol);
    validate<"[!#1]~[#7](~[!#1])~[!#1]">(mol);
}

template void RDMACCS_part_24<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_24<RDKit::ROMol>(RDKit::ROMol &mol);
