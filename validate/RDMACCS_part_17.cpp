#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_17(Mol &mol)
{
    // SMARTS 81 - 85
    validate<"[!#6;!#1]~1~*~*~*~*~1">(mol);
    validate<"[NH2]">(mol);
    validate<"[#6]~[#7](~[#6])~[#6]">(mol);
    validate<"[C;H2,H3][!#6;!#1][C;H2,H3]">(mol);
    validate<"[F,Cl,Br,I]!@*@*">(mol);
}

template void RDMACCS_part_17<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_17<RDKit::ROMol>(RDKit::ROMol &mol);
