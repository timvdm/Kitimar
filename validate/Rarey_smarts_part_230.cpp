#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_230(Mol &mol)
{
    // SMARTS 1146 - 1150
    validate<"c:1:c:c:c:c(:c:1-[#7](-[#1])-[!$([#6]=[#8])])-[#6](-[#6]:[#6])=[#8]">(mol);
    validate<"c:1:c:c:c:c:2:c:1:c:c:3:c(:n:2):n:c:4:c(:c:3-[#7]):c:c:c:c:4">(mol);
    validate<"c:1:c:c:c:c:c:1-[#6](=[#8])-[#7](-[#1])-[#7]=[#6]-3-c:2:c:c:c:c:c:2-c:4:c:c:c:c:c-3:4">(mol);
    validate<"c:1:c:c:c:c:c:1-[#7](-[#1])-[#6](=[#16])-[#7](-[#1])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-c:2:c:c:c:c:c:2">(mol);
    validate<"c:1:c:c:c:c:c:1-[#7](-[#1])-[#6](=[#16])-[#7](-[#1])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-c:2:c:c:c:c:c:2">(mol);
}

template void Rarey_smarts_part_230<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_230<RDKit::ROMol>(RDKit::ROMol &mol);
