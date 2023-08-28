#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_215(Mol &mol)
{
    // SMARTS 1071 - 1075
    validate<"c:1-2:c(:c:c:c:c:1)-[#6](=[#8])-[#6](=[#6])-[#6]-2=[#8]">(mol);
    validate<"c:1-2:c(:c:c:c:c:1-[#6](-[#1])(-[#1])-[#6](-[#1])=[#6](-[#1])-[#1])-[#6](=[#6](-[#6](=[#8])-[#7](-[#1])-[#6]:[#6])-[#6](=[#8])-[#8]-2)-[#1]">(mol);
    validate<"c:1-2:c:c-3:c(:c:c:1-[#8]-[#6]-[#8]-2)-[#6]-[#6]-3">(mol);
    validate<"c:1-3:c(:c(:c(:c(:c:1)-[#8]-[#6]-[#1])-[#1])-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])-[#1])-[#6](=[#8])-[#8]-3">(mol);
    validate<"c:1-3:c(:c(:c(:c(:c:1-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])-[#1])-[#6](=[#7]-[#7](-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#6](=[#8])-[#8]-[#1])-[#1])-[#1])-c:4:c-3:c(:c(:c(:c:4-[#1])-[#8]-[#6](-[#1])-[#1])-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_215<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_215<RDKit::ROMol>(RDKit::ROMol &mol);
