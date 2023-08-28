#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_10(Mol &mol)
{
    // SMARTS 46 - 46
    validate<"[#6]=,:1[#6]=,:[#6][*D2][#6]=,:1">(mol);
}

template void Agrafiotis_ABCD_part_10<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_10<RDKit::ROMol>(RDKit::ROMol &mol);
