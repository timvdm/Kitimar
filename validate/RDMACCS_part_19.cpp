#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_19(Mol &mol)
{
    // SMARTS 91 - 95
    validate<"[!#6;!#1]~[CH3]">(mol);
    validate<"[!#6;!#1]~[#7]">(mol);
    validate<"[#7]~*~*~[#8]">(mol);
    validate<"*~1~*~*~*~*~1">(mol);
    validate<"[#7]~*~*~*~[#8]">(mol);
}

template void RDMACCS_part_19<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_19<RDKit::ROMol>(RDKit::ROMol &mol);
