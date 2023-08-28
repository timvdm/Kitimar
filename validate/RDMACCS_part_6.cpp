#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_6(Mol &mol)
{
    // SMARTS 26 - 30
    validate<"[I]">(mol);
    validate<"[!#6;!#1]~[CH2]~[!#6;!#1]">(mol);
    validate<"[#15]">(mol);
    validate<"[#6]~[!#6;!#1](~[#6])(~[#6])~[!#1]">(mol);
    validate<"[!#6;!#1]~[F,Cl,Br,I]">(mol);
}

template void RDMACCS_part_6<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_6<RDKit::ROMol>(RDKit::ROMol &mol);
