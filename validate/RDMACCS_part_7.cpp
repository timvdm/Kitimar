#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_7(Mol &mol)
{
    // SMARTS 31 - 35
    validate<"[#6]~[#16]~[#7]">(mol);
    validate<"[#7]~[#16]">(mol);
    validate<"[CH2]=*">(mol);
    validate<"[Li,Na,K,Rb,Cs,Fr]">(mol);
    validate<"[#16R]">(mol);
}

template void RDMACCS_part_7<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_7<RDKit::ROMol>(RDKit::ROMol &mol);
