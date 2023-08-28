#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_4(Mol &mol)
{
    // SMARTS 16 - 20
    validate<"[#6]#[#6]">(mol);
    validate<"[#5,Al,Ga,In,Tl]">(mol);
    validate<"*~1~*~*~*~*~*~*~1">(mol);
    validate<"[#14]">(mol);
    validate<"[#6]=[#6](~[!#6;!#1])~[!#6;!#1]">(mol);
}

template void RDMACCS_part_4<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_4<RDKit::ROMol>(RDKit::ROMol &mol);
