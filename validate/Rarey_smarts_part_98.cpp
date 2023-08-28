#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_98(Mol &mol)
{
    // SMARTS 486 - 490
    validate<"[#7]-2(-c:1:c:c:c:c:c:1)-[#7]=[#6](-[#7](-[#1])-[#6]=[#8])-[#6](-[#1])(-[#1])-[#6]-2=[#8]">(mol);
    validate<"[#7]-2(-c:1:c:c:c:c:c:1-[#6](-[#1])-[#1])-[#6](=[#16])-[#7](-[#6](-[#1])(-[#1])-[!#1]:[!#1]:[!#1]:[!#1]:[!#1])-[#6](-[#1])(-[#1])-[#6]-2=[#8]">(mol);
    validate<"[#7]-2-[#16]-[#6]-1=[#6](-[#6]:[#6]-[#7]-[#6]-1)-[#6]-2=[#16]">(mol);
    validate<"[#7]-2-[#6]=[#6](-[#6]=[#8])-[#6](-c:1:c:c:c(:c:c:1)-[#7](-[#6](-[#1])-[#1])-[#6](-[#1])-[#1])-[#6]~3=[#6]-2~[#7]~[#6](~[#16])~[#7]~[#6]~3~[#7]">(mol);
    validate<"[#7]-2-c:1:c:c:c:c:c:1-[#6](=[#7])-c:3:c-2:c:c:c:c:3">(mol);
}

template void Rarey_smarts_part_98<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_98<RDKit::ROMol>(RDKit::ROMol &mol);
