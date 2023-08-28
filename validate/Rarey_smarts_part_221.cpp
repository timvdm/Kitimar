#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_221(Mol &mol)
{
    // SMARTS 1101 - 1105
    validate<"c:1:c(:c:c:c:c:1)-[#6](-[#1])(-[#1])-[#7](-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#8]-[#1])-[#1])-[#1]">(mol);
    validate<"c:1:c(:c:c:c:c:1)-[#6](-[#1])-[#7]-[#6](=[#8])-[#6](-[#7](-[#1])-[#6](-[#1])-[#1])=[#6](-[#1])-[#6](=[#8])-c:2:c:c:c(:c:c:2)-[#8]-[#6](-[#1])-[#1]">(mol);
    validate<"c:1:c(:c:c:c:c:1)-[#6](=[#8])-[#6](-[#1])=[#6]-3-[#6](=[#8])-[#7](-[#1])-[#6](=[#8])-[#6](=[#6](-[#1])-c:2:c:c:c:c:c:2)-[#7]-3-[#1]">(mol);
    validate<"c:1:c(:c:c:c:c:1)-[#6](=[#8])-[#7](-[#1])-c:2:c(:c:c:c:c:2)-[#6](=[#8])-[#7](-[#1])-[#7](-[#1])-c:3:n:c:c:s:3">(mol);
    validate<"c:1:c(:c:c:c:c:1)-[#6]-4=[#7]-[#7]:2:[#6](:[#7+]:c:3:c:2:c:c:c:c:3)-[#16]-[#6;X4]-4">(mol);
}

template void Rarey_smarts_part_221<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_221<RDKit::ROMol>(RDKit::ROMol &mol);
