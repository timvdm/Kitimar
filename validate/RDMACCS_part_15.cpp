#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_15(Mol &mol)
{
    // SMARTS 71 - 75
    validate<"[#16]=*">(mol);
    validate<"[CH3]~*~[CH3]">(mol);
    validate<"[!#1]!@[#7]@[!#1]">(mol);
    validate<"[#6]=[#6](~[!#1])~[!#1]">(mol);
    validate<"[#7]~*~[#7]">(mol);
}

template void RDMACCS_part_15<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_15<RDKit::ROMol>(RDKit::ROMol &mol);
