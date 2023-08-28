#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_28(Mol &mol)
{
    // SMARTS 136 - 140
    validate<"[O;!H0]">(mol);
    validate<"[#8]">(mol);
    validate<"[CH3]">(mol);
    validate<"[#7]">(mol);
    validate<"*@*!@[#8]">(mol);
}

template void RDMACCS_part_28<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_28<RDKit::ROMol>(RDKit::ROMol &mol);
