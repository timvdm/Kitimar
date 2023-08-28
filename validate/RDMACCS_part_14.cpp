#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_14(Mol &mol)
{
    // SMARTS 66 - 70
    validate<"[!#6;!#1;!H0]~[!#6;!#1;!H0]">(mol);
    validate<"[!#6;!#1]~[!#6;!#1;!H0]">(mol);
    validate<"[!#6;!#1]~[#7]~[!#6;!#1]">(mol);
    validate<"[#7]~[#8]">(mol);
    validate<"[#8]~*~*~[#8]">(mol);
}

template void RDMACCS_part_14<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_14<RDKit::ROMol>(RDKit::ROMol &mol);
