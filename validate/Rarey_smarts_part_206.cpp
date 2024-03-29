#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_206(Mol &mol)
{
    // SMARTS 1026 - 1030
    validate<"c:1(:c(:c(:c(:c(:c:1-[$([#1]),$([#6](-[#1])-[#1])])-[#1])-[#8]-[#6](-[#1])-[#1])-[#6](-[#1])(-[#1])-[$([#7](-[#1])-[#6](=[#8])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#6](-[#1])-[#1]),$([#6](-[#1])(-[#6](-[#1])-[#1])-[#7](-[#1])-[#6](=[#16])-[#7]-[#1])])-[#1])-[#8]-[#6](-[#1])-[#1]">(mol);
    validate<"c:1(:c(:c(:c(:o:1)-[#6](-[#1])-[#1])-[#1])-[#1])-[#6](-[#1])(-[#8]-[#1])-[#6]#[#6]-[#6;X4]">(mol);
    validate<"c:1(:c(:c(:c(:o:1)-[#6](-[#1])-[#1])-[#6](-[#1])(-[#1])-[#8]-[#6]:[#6])-[#1])-[#6](=[#8])-[#8]-[#1]">(mol);
    validate<"c:1(:c(:c(:c(:o:1)-[$([#1]),$([#6](-[#1])-[#1])])-[#1])-[#1])-[#6](-[$([#1]),$([#6](-[#1])-[#1])])=[#7]-[#7](-[#1])-c:2:c:c:n:c:c:2">(mol);
    validate<"c:1(:c(:c(:c(:o:1)-[$([#1]),$([#6](-[#1])-[#1])])-[#1])-[#1])-[#6](-[$([#1]),$([#6](-[#1])-[#1])])=[#7]-[#7](-[#1])-c:2:n:c:c:s:2">(mol);
}

template void Rarey_smarts_part_206<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_206<RDKit::ROMol>(RDKit::ROMol &mol);
