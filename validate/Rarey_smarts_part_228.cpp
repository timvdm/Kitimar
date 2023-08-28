#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_228(Mol &mol)
{
    // SMARTS 1136 - 1140
    validate<"c:1:c:c-2:c(:c:c:1)-[#6]-[#16]-c3c(-[#6]-2=[#6])ccs3">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#6]-[#6](-c:3:c(-[#16]-2):c(:c(-[#1]):c(:c:3-[#1])-[$([#1]),$([#8]),$([#16;X2]),$([#6;X4]),$([#7](-[$([#1]),$([#6;X4])])-[$([#1]),$([#6;X4])])])-[#1])-[#7](-[$([#1]),$([#6;X4])])-[$([#1]),$([#6;X4])]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#6]=[#6]-[#6](-[#7]-2-[#6](=[#8])-[#7](-[#1])-c:3:c:c(:c(:c:c:3)-[#8]-[#6](-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])(-[#6](-[#1])-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#7](-[#6](-[#8]-[#6]-2)(-[#6](=[#8])-[#8]-[#1])-[#6](-[#1])-[#1])-[#6](=[#8])-[#6](-[#1])-[#1]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#7]=[#6]-[#6]-2=[#7;!R]">(mol);
}

template void Rarey_smarts_part_228<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_228<RDKit::ROMol>(RDKit::ROMol &mol);
