#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_217(Mol &mol)
{
    // SMARTS 1081 - 1085
    validate<"c:12:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1])c(c(-[#6]:[#6])n2-!@[#6]:[#6])-[#6](-[#1])-[#1]">(mol);
    validate<"c:12:c(:c:c:c:n:1)c(c(-[#6](=[#8])~[#8;X1])s2)-[#7](-[#1])-[#1]">(mol);
    validate<"c:1:2(:c(:c(:c(:o:1)-[#6])-[#1])-[#1])-[#6](=[#8])-[#7](-[#1])-[#6]:[#6](-[#1]):[#6](-[#1]):[#6](-[#1]):[#6](-[#1]):[#6]:2-[#6](=[#8])-[#8]-[#1]">(mol);
    validate<"c:1:2:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1]):[!#6&!#1]:[#6](:[#6]:2-[#6](-[#1])=[#7]-[#7](-[#1])-[$([#6]:1:[#7]:[#6]:[#6](-[#1]):[#16]:1),$([#6]:[#6](-[#1]):[#6]-[#1]),$([#6]:[#7]:[#6]:[#7]:[#6]:[#7]),$([#6]:[#7]:[#7]:[#7]:[#7])])-[$([#1]),$([#8]-[#1]),$([#6](-[#1])-[#1])]">(mol);
    validate<"c:1:2:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1]):o:c:3:c(-[#1]):c(:c(-[#8]-[#6](-[#1])-[#1]):c(:c:2:3)-[#1])-[#7](-[#1])-[#6](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_217<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_217<RDKit::ROMol>(RDKit::ROMol &mol);