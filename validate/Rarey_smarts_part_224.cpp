#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_224(Mol &mol)
{
    // SMARTS 1116 - 1120
    validate<"c:1:c-3:c(:c:c:c:1)-[#6]:2:[#7]:[!#1]:[#6]:[#6]:[#6]:2-[#6]-3=[#8]">(mol);
    validate<"c:1:c-3:c(:c:c:c:1)-[#7](-c:2:c:c:c:c:c:2-[#8]-3)-[#6](-[#1])(-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"c:1:c:2:c(:c:c:c:1):c(:c:3:c(:c:2):c:c:c:c:3)-[#6]=[#7]-[#7](-[#1])-c:4:c:c:c:c:c:4">(mol);
    validate<"c:1:c:2:c(:c:c:c:1):n:c:3:c(:c:2-[#7]):c:c:c:c:3">(mol);
    validate<"c:1:c:4:c(:c:c2:c:1nc(n2-[#1])-[#6]-[#8]-[#6](=[#8])-c:3:c:c(:c:c(:c:3)-[#7](-[#1])-[#1])-[#7](-[#1])-[#1]):c:c:c:c:4">(mol);
}

template void Rarey_smarts_part_224<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_224<RDKit::ROMol>(RDKit::ROMol &mol);
