#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_96(Mol &mol)
{
    // SMARTS 476 - 480
    validate<"[#7]-1=[#6]-[#6](-[#6](-[#7]-1)=[#16])=[#6]">(mol);
    validate<"[#7]-2(-[#6](-[#1])-[#1])-[#6](=[#16])-[#7](-[#1])-[#6](=[#6](-[#1])-c:1:c:c:c:c(:c:1)-[Br])-[#6]-2=[#8]">(mol);
    validate<"[#7]-2(-[#6](=[#8])-c:1:c(:c(:c(:c(:c:1-[#1])-[#6](=[#8])-[#8]-[#1])-[#1])-[#1])-[#6]-2=[#8])-c:3:c(:c:c(:c(:c:3)-[#1])-[#8])-[#1]">(mol);
    validate<"[#7]-2(-c:1:c:c:c(:c:c:1)-[#8]-[#6](-[#1])-[#1])-[#6](=[#8])-[#6](=[#6]-[#6](=[#7]-2)-n:3:c:n:c:c:3)-[#6]#[#7]">(mol);
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#6](=[#7]-[#6]=[#8])-[#16]-[#6](-[#1])(-[#1])-[#6]-2=[#8]">(mol);
}

template void Rarey_smarts_part_96<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_96<RDKit::ROMol>(RDKit::ROMol &mol);
