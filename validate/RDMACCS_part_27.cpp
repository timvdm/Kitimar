#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_27(Mol &mol)
{
    // SMARTS 131 - 135
    validate<"[F,Cl,Br,I]">(mol);
    validate<"[#7]!:*:*">(mol);
    validate<"[#8]=*">(mol);
    validate<"[!#6;R]">(mol);
    validate<"[!#6;!#1]~[CH2]~[!#1]">(mol);
}

template void RDMACCS_part_27<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_27<RDKit::ROMol>(RDKit::ROMol &mol);
