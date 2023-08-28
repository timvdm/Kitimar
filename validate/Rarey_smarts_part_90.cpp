#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_90(Mol &mol)
{
    // SMARTS 446 - 450
    validate<"[#7](-[#1])(-[#1])-c:1:c(:c(:c(:n:c:1-[#1])-[#8]-c:2:c:c:c:c:c:2)-[#1])-[#1]">(mol);
    validate<"[#7](-[#1])(-[#1])-c:1:c(:c(:c(:s:1)-[!#1])-[!#1])-[#6]=[#8]">(mol);
    validate<"[#7](-[#1])(-[#1])-c:1:c(:c(:c(:s:1)-[#7](-[#1])-[#6](=[#8])-c:2:c:c:c:c:c:2)-[#6]#[#7])-[#6]:3:[!#1]:[!#1]:[!#1]:[!#1]:[!#1]:3">(mol);
    validate<"[#7](-[#1])(-[#1])-c:1:c(:c:c:c:n:1)-[#8]-[#6](-[#1])(-[#1])-[#6]:[#6]">(mol);
    validate<"[#7](-[#1])(-[#1])-c:1:c:c:c(:c:c:1-[#8]-[#1])-[#16](=[#8])(=[#8])-[#8]-[#1]">(mol);
}

template void Rarey_smarts_part_90<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_90<RDKit::ROMol>(RDKit::ROMol &mol);
