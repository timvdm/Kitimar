#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_97(Mol &mol)
{
    // SMARTS 481 - 485
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#6](=[#8])-[#16]-[#6](-[#1])(-[#1])-[#6]-2=[#16]">(mol);
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#6](=[#8])-[#16]-[#6](-[#1])(-[#6](-[#1])(-[#1])-[#6](=[#8])-[#7](-[#1])-[#6]:[#6])-[#6]-2=[#8]">(mol);
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#6](=[#8])-[#6](=[#6]-[#6](=[#7]-2)-[#6]#[#7])-[#6]#[#7]">(mol);
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#7]=[#6](-[#6](-[#1])-[#1])-[#6](-[#1])(-[#16]-[#6])-[#6]-2=[#8]">(mol);
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#7]=[#6](-[#6]=[#8])-[#6;X4]-[#6]-2=[#8]">(mol);
}

template void Rarey_smarts_part_97<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_97<RDKit::ROMol>(RDKit::ROMol &mol);
