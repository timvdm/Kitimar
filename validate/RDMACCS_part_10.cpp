#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_10(Mol &mol)
{
    // SMARTS 46 - 50
    validate<"[#8]~[!#6;!#1](~[#8])(~[#8])">(mol);
    validate<"[!+0]">(mol);
    validate<"[#6]=[#6](~[#6])~[#6]">(mol);
    validate<"[#6]~[#16]~[#8]">(mol);
    validate<"[#7]~[#7]">(mol);
}

template void RDMACCS_part_10<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_10<RDKit::ROMol>(RDKit::ROMol &mol);
