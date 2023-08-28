#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_107(Mol &mol)
{
    // SMARTS 531 - 535
    validate<"[#8]=[#6]-1-[#7]-[#7]-[#6](=[#7]-[#6]-1=[#6]-[#1])-[!#1]:[!#1]">(mol);
    validate<"[#8]=[#6]-2-[#16]-c:1:c(:c(:c:c:c:1)-[#8]-[#6](-[#1])-[#1])-[#8]-2">(mol);
    validate<"[#8]=[#6]-2-[#6](=!@[#7]-[#7])-c:1:c:c:c:c:c:1-[#7]-2">(mol);
    validate<"[#8]=[#6]-3-[#6](=!@[#6](-[#1])-c:1:c:n:c:c:1)-c:2:c:c:c:c:c:2-[#7]-3">(mol);
    validate<"[#8]=[#6]-3-c:1:c(:c:c:c:c:1)-[#6]-2=[#6](-[#8]-[#1])-[#6](=[#8])-[#7]-c:4:c-2:c-3:c:c:c:4">(mol);
}

template void Rarey_smarts_part_107<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_107<RDKit::ROMol>(RDKit::ROMol &mol);
