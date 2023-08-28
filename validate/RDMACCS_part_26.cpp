#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_26(Mol &mol)
{
    // SMARTS 126 - 130
    validate<"[$([!#1]~[CH2]~*~*~[CH2]~[!#1]),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$([!#1]~[CH2]~[R]1@[R]@[CH2;R]1)]">(mol);
    validate<"[!#6;!#1]~[!#6;!#1]">(mol);
    validate<"[!#6;!#1;!H0]">(mol);
    validate<"[#8]~*~[CH2]~[!#1]">(mol);
    validate<"*@*!@[#7]">(mol);
}

template void RDMACCS_part_26<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_26<RDKit::ROMol>(RDKit::ROMol &mol);
