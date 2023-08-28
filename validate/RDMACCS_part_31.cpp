#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_31(Mol &mol)
{
    // SMARTS 151 - 155
    validate<"[#6]=[#8]">(mol);
    validate<"[!#1]!@[CH2]!@[!#1]">(mol);
    validate<"[#7]~[!#1](~[!#1])~[!#1]">(mol);
    validate<"[#6]-[#8]">(mol);
    validate<"[#6]-[#7]">(mol);
}

template void RDMACCS_part_31<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_31<RDKit::ROMol>(RDKit::ROMol &mol);
