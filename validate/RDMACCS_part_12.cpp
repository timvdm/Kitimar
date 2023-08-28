#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_12(Mol &mol)
{
    // SMARTS 56 - 60
    validate<"[!#6;!#1]~[#16]~[!#6;!#1]">(mol);
    validate<"[#16]!:*:*">(mol);
    validate<"[#16]=[#8]">(mol);
    validate<"[!#1]~[#16](~[!#1])~[!#1]">(mol);
    validate<"*@*!@*@*">(mol);
}

template void RDMACCS_part_12<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_12<RDKit::ROMol>(RDKit::ROMol &mol);
