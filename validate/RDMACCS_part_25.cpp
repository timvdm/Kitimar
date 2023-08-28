#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_25(Mol &mol)
{
    // SMARTS 121 - 125
    validate<"[#8]~[#6]~[#8]">(mol);
    validate<"[!#6;!#1]~[!#6;!#1]">(mol);
    validate<"[!#1]!@[#8]!@[!#1]">(mol);
    validate<"*@*!@[#8]">(mol);
    validate<"[$([!#1]~[CH2]~*~*~*~[CH2]~[!#1]),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$([!#1]~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$([!#1]~[CH2]~*~[R]1@[R]@[CH2;R]1)]">(mol);
}

template void RDMACCS_part_25<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_25<RDKit::ROMol>(RDKit::ROMol &mol);
