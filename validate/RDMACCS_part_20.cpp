#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_20(Mol &mol)
{
    // SMARTS 96 - 100
    validate<"[!#6;!#1]~1~*~*~*~*~*~1">(mol);
    validate<"[#6]=[#6]">(mol);
    validate<"[!#1]~[CH2]~[#7]">(mol);
    validate<"[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]">(mol);
    validate<"[!#6;!#1]~[#8]">(mol);
}

template void RDMACCS_part_20<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_20<RDKit::ROMol>(RDKit::ROMol &mol);
