#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_8(Mol &mol)
{
    // SMARTS 36 - 40
    validate<"[#7]~[#6](~[#8])~[#7]">(mol);
    validate<"[#7]~[#6](~[#6])~[#7]">(mol);
    validate<"[#8]~[#16](~[#8])~[#8]">(mol);
    validate<"[#16]-[#8]">(mol);
    validate<"[#6]#[#7]">(mol);
}

template void RDMACCS_part_8<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_8<RDKit::ROMol>(RDKit::ROMol &mol);
