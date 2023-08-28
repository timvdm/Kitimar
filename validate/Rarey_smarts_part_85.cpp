#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_85(Mol &mol)
{
    // SMARTS 421 - 425
    validate<"[#6][F,Cl,Br,I]">(mol);
    validate<"[#6][OX2H]">(mol);
    validate<"[#6]~!@[#8]">(mol);
    validate<"[#6]~1~3~[#7](-[#6]:[#6])~[#6]~[#6]~[#6]~[#6]~1~[#6]~2~[#7]~[#6]~[#6]~[#6]~[#7+]~2~[#7]~3">(mol);
    validate<"[#6]~1~[#6](~[#7]~[#7]~[#6](~[#6](-[#1])-[#1])~[#6](-[#1])-[#1])~[#7]~[#16]~[#6]~1">(mol);
}

template void Rarey_smarts_part_85<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_85<RDKit::ROMol>(RDKit::ROMol &mol);
