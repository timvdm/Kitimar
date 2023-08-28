#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_13(Mol &mol)
{
    // SMARTS 61 - 65
    validate<"[#7]=[#8]">(mol);
    validate<"*@*!@[#16]">(mol);
    validate<"c:n">(mol);
    validate<"[#6]~[#6](~[#6])(~[#6])~[!#1]">(mol);
    validate<"[!#6;!#1]~[#16]">(mol);
}

template void RDMACCS_part_13<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_13<RDKit::ROMol>(RDKit::ROMol &mol);
