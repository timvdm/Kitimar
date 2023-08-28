#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_5(Mol &mol)
{
    // SMARTS 21 - 25
    validate<"*~1~*~*~1">(mol);
    validate<"[#7]~[#6](~[#8])~[#8]">(mol);
    validate<"[#7]-[#8]">(mol);
    validate<"[#7]~[#6](~[#7])~[#7]">(mol);
    validate<"[#6]=;@[#6](@*)@*">(mol);
}

template void RDMACCS_part_5<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_5<RDKit::ROMol>(RDKit::ROMol &mol);
