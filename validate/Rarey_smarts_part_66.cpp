#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_66(Mol &mol)
{
    // SMARTS 326 - 330
    validate<"[#6]-1(-[#6](-[#1])(-[#1])-[#6]-1(-[#1])-[#1])(-[#6](=[#8])-[#7](-[#1])-c:2:c:c:c(:c:c:2)-[#8]-[#6](-[#1])(-[#1])-[#8])-[#16](=[#8])(=[#8])-[#6]:[#6]">(mol);
    validate<"[#6]-1(-[#6](-[#6]=[#6]-[!#6&!#1]-1)=[#6])=[!#6&!#1]">(mol);
    validate<"[#6]-1(-[#6](=[#6](-[#6]#[#7])-[#6](~[#8])~[#7]~[#6]-1~[#8])-[#6](-[#1])-[#1])=[#6](-[#1])-[#6]:[#6]">(mol);
    validate<"[#6]-1(-[#6](=[#6]-[#6]=[#6]-[#6]=[#6]-1)-[#7]-[#1])=[#7]-[#6]">(mol);
    validate<"[#6]-1(-[#6](=[#8])-[#6](-[#1])(-[#1])-[#6]-[#6](-[#1])(-[#1])-[#6]-1=[#8])=[#6](-[#7]-[#1])-[#6]=[#8]">(mol);
}

template void Rarey_smarts_part_66<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_66<RDKit::ROMol>(RDKit::ROMol &mol);
