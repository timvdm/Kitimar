#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_9(Mol &mol)
{
    // SMARTS 41 - 45
    validate<"F">(mol);
    validate<"[!C;!c;!#1;!H0]~*~[!C;!c;!#1;!H0]">(mol);
    validate<"[#6]=[#6]~[#7]">(mol);
    validate<"Br">(mol);
    validate<"[#16]~*~[#7]">(mol);
}

template void RDMACCS_part_9<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_9<RDKit::ROMol>(RDKit::ROMol &mol);
