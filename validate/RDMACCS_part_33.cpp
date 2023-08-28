#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_33(Mol &mol)
{
    // SMARTS 161 - 162
    validate<"[#8]">(mol);
    validate<"[R]">(mol);
}

template void RDMACCS_part_33<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_33<RDKit::ROMol>(RDKit::ROMol &mol);
