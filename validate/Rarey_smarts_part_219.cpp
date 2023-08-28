#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_219(Mol &mol)
{
    // SMARTS 1091 - 1095
    validate<"c:1:2:c(:c(:c(:c(:c:1:c(:c(:c(:c:2-[#1])-[#8]-[#1])-[#6](=[#8])-[#7](-[#1])-[#7]=[#6])-[#1])-[#1])-[#1])-[#1])-[#1]">(mol);
    validate<"c:1:2:c:3:c(:c(-[#8]-[#1]):c(:c:1:c(:c:n:2-[#6])-[#6]=[#8])-[#1]):n:c:n:3">(mol);
    validate<"c:1:2:c:c:c:c(:c:1:c(:c:c:c:2)-[$([#8]-[#1]),$([#7](-[#1])-[#1])])-[#6](-[#6])=[#8]">(mol);
    validate<"c:1:2:n:c(:c(:n:c:1:[#6]:[#6]:[#6]:[!#1]:2)-[#6](-[#1])=[#6](-[#8]-[#1])-[#6])-[#6](-[#1])=[#6](-[#8]-[#1])-[#6]">(mol);
    validate<"c:1:3:c(:c(:c(:c(:c:1-[#1])-[#1])-[#7](-[#1])-[#6](-[#1])(-[#1])-c:2:c:c:c:c:c:2)-[#1]):n:c(-[#1]):n:3-[#6]">(mol);
}

template void Rarey_smarts_part_219<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_219<RDKit::ROMol>(RDKit::ROMol &mol);
