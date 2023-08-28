#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_103(Mol &mol)
{
    // SMARTS 511 - 515
    validate<"[#7]~[D3]">(mol);
    validate<"[#8,#16]">(mol);
    validate<"[#8](-[#1])-[#6](=[#8])-c:1:c(:c(:c(:c(:c:1-[#8]-[#1])-[#1])-c:2:c(-[#1]):c(:c(:o:2)-[#6](-[#1])=[#6](-[#6]#[#7])-c:3:n:c:c:n:3)-[#1])-[#1])-[#1]">(mol);
    validate<"[#8](-[#1])-[#6](=[#8])-c:1:c:c(:c:c:c:1)-[#6]:[!#1]:[#6]-[#6](-[#1])=[#6]-2-[#6](=[!#6&!#1])-[#7]-[#6](=[!#6&!#1])-[!#6&!#1]-2">(mol);
    validate<"[#8](-[#1])-[#6](=[#8])-c:1:c:c:c(:c:c:1)-[#7]-[#7]=[#6](-[#1])-[#6]:2:[#6](:[#6](:[#6](:[!#1]:2)-c:3:c:c:c:c:c:3)-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_103<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_103<RDKit::ROMol>(RDKit::ROMol &mol);
