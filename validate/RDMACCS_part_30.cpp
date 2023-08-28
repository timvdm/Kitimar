#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_30(Mol &mol)
{
    // SMARTS 146 - 150
    validate<"[C;H3,H4]">(mol);
    validate<"[!#1]!@*@*!@[!#1]">(mol);
    validate<"[#7;!H0]">(mol);
    validate<"[#8]~[#6](~[#6])~[#6]">(mol);
    validate<"[!#6;!#1]~[CH2]~[!#1]">(mol);
}

template void RDMACCS_part_30<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_30<RDKit::ROMol>(RDKit::ROMol &mol);
