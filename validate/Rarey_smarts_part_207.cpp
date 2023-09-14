#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_207(Mol &mol)
{
    // SMARTS 1031 - 1035
    validate<"c:1(:c(:c(:c(:o:1)-[$([#1]),$([#6](-[#1])-[#1])])-[#1])-[#1])-[#6](=[#8])-[#7](-[#1])-[#7]=[#6](-[$([#1]),$([#6](-[#1])-[#1])])-c:2:c:c:c:c(:c:2)-[*]-[*]-[*]-c:3:c:c:c:o:3">(mol);
    validate<"c:1(:c(:c(:c(:s:1)-[#1])-[#1])-[$([#1]),$([#6](-[#1])-[#1])])-[#6](-[#1])=[#7]-[#7](-[#1])-c:2:c:c:c:c:c:2">(mol);
    validate<"c:1(:c(:c(:c(:s:1)-[$([#1]),$([#6](-[#1])-[#1])])-[#1])-[#1])-[#6](-[$([#1]),$([#6](-[#1])-[#1])])-[#6](=[#8])-[#7](-[#1])-c:2:n:c:c:s:2">(mol);
    validate<"c:1(:c(:c-2:c(:c(:c:1-[#1])-[#1])-[#7](-[#6](-[#7]-2-[#1])=[#8])-[#1])-[#1])-[#7](-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"c:1(:c(:c-2:c(:c(:c:1-[#8]-[#6](-[#1])-[#1])-[#1])-[#6]=[#6]-[#6](-[#1])-[#16]-2)-[#1])-[#8]-[#6](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_207<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_207<RDKit::ROMol>(RDKit::ROMol &mol);