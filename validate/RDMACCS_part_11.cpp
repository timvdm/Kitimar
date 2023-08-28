#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_11(Mol &mol)
{
    // SMARTS 51 - 55
    validate<"[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]">(mol);
    validate<"[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]">(mol);
    validate<"[#8]~[#16]~[#8]">(mol);
    validate<"[#8]~[#7](~[#8])~[#6]">(mol);
    validate<"[#8R]">(mol);
}

template void RDMACCS_part_11<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_11<RDKit::ROMol>(RDKit::ROMol &mol);
