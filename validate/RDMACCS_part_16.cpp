#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_16(Mol &mol)
{
    // SMARTS 76 - 80
    validate<"[#6]=[#7]">(mol);
    validate<"[#7]~*~*~[#7]">(mol);
    validate<"[#7]~*~*~*~[#7]">(mol);
    validate<"[#16]~*(~[!#1])~[!#1]">(mol);
    validate<"[!#1]~[CH2]~[!#6;!#1;!H0]">(mol);
}

template void RDMACCS_part_16<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_16<RDKit::ROMol>(RDKit::ROMol &mol);
