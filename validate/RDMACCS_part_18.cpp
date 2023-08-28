#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_18(Mol &mol)
{
    // SMARTS 86 - 90
    validate<"[#16]">(mol);
    validate<"[#8]~*~*~*~[#8]">(mol);
    validate<"[$([!#6;!#1;!H0]~*~*~[CH2]~[!#1]),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]">(mol);
    validate<"[$([!#6;!#1;!H0]~*~*~*~[CH2]~[!#1]),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]">(mol);
    validate<"[#8]~[#6](~[#7])~[#6]">(mol);
}

template void RDMACCS_part_18<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_18<RDKit::ROMol>(RDKit::ROMol &mol);
