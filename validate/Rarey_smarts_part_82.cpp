#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_82(Mol &mol)
{
    // SMARTS 406 - 410
    validate<"[#6]:[#7$(a1aaaaa1)]">(mol);
    validate<"[#6]:[#8$(a1aaaa1)]">(mol);
    validate<"[#6]:[#8$(a1aaaaa1)]">(mol);
    validate<"[#6]=!@[#6](-[!#1])-@[#6](=!@[!#6&!#1])-@[#6](=!@[#6])-[!#1]">(mol);
    validate<"[#6]=[#6](-[#6]#[#7])-[#6](=[#7]-[#1])-[#7]-[#7]">(mol);
}

template void Rarey_smarts_part_82<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_82<RDKit::ROMol>(RDKit::ROMol &mol);
