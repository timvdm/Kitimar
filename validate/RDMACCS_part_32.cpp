#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_32(Mol &mol)
{
    // SMARTS 156 - 160
    validate<"[#8]">(mol);
    validate<"[C;H3,H4]">(mol);
    validate<"[#7]">(mol);
    validate<"a">(mol);
    validate<"*~1~*~*~*~*~*~1">(mol);
}

template void RDMACCS_part_32<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_32<RDKit::ROMol>(RDKit::ROMol &mol);
