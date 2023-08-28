#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_104(Mol &mol)
{
    // SMARTS 516 - 520
    validate<"[#8](-[#1])-[#6](=[#8])-c:1:c:c:c:c(:c:1)-[#6]:[!#1]:[#6]-[#6]=[#7]-[#7](-[#1])-[#6](=[#8])-[#6](-[#1])(-[#1])-[#8]">(mol);
    validate<"[#8](-[#1])-[#6]:1:[#6](:[#6]:[!#1]:[#6](:[#7]:1)-[#7](-[#1])-[#1])-[#6](-[#1])(-[#1])-[#6](=[#8])-[#8]">(mol);
    validate<"[#8](-[#1])-c:1:n:c(:c:c:c:1)-[#8]-[#1]">(mol);
    validate<"[#8](-c:1:c:c:c:c:c:1)-c:3:c:c:2:n:o:n:c:2:c:c:3">(mol);
    validate<"[#8]-1-[#6](-[#16]-c:2:c-1:c:c:c(:c:2)-[$([#7]),$([#8])])=[$([#8]),$([#16])]">(mol);
}

template void Rarey_smarts_part_104<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_104<RDKit::ROMol>(RDKit::ROMol &mol);
