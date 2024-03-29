#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_58(Mol &mol)
{
    // SMARTS 286 - 290
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-[#6]-2=[#6]-c:1:c(:c:c:c:c:1)-[#6]-2(-[#1])-[#1]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-c:1:c(-[#1]):c(:c(:o:1)-[#6](-[#1])=[#6]-[#6]#[#7])-[#1]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#6](-[#1])=[#7]-[#7]-[$([#6](=[#8])-[#6](-[#1])(-[#1])-[#16]-[#6]:[#7]),$([#6](=[#8])-[#6](-[#1])(-[#1])-[!#1]:[!#1]:[#7]),$([#6](=[#8])-[#6]:[#6]-[#8]-[#1]),$([#6]:[#7]),$([#6](-[#1])(-[#1])-[#6](-[#1])-[#8]-[#1])])-[#1])-[#1]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#6](-[#1])=[#7]-[#7]=[#6](-[#6])-[#6]:[#6])-[#1])-[#1]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-c:1:c:c(:c(:c(:c:1)-[$([#1]),$([#6](-[#1])-[#1]),$([#8]-[#6](-[#1])(-[#1])-[#6](-[#1])-[#1])])-[#7])-[#1]">(mol);
}

template void Rarey_smarts_part_58<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_58<RDKit::ROMol>(RDKit::ROMol &mol);
