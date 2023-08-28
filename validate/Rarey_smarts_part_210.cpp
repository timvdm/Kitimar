#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_210(Mol &mol)
{
    // SMARTS 1046 - 1050
    validate<"c:1(:c:c:c(:c:c:1)-[#6](-[#1])-[#1])-c:2:c(:s:c(:n:2)-[#7](-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1]">(mol);
    validate<"c:1(:c:c:c(:c:c:1)-[#6]-3=[#6]-[#6](-c2cocc2-[#6](=[#6]-3)-[#8]-[#1])=[#8])-[#16]-[#6](-[#1])-[#1]">(mol);
    validate<"c:1(:c:c:c(:c:c:1)-[#6]=[#7]-[#7])-[#8]-[#1]">(mol);
    validate<"c:1(:c:c:c:c:c:1)-[#6](-[#1])=!@[#6]-3-[#6](=[#8])-c:2:c:c:c:c:c:2-[#16]-3">(mol);
    validate<"c:1(:c:c:c:c:c:1)-[#7](-[#1])-[#6](=[#16])-[#7](-[#1])-[#7]=[#6]-c:2:c:n:c:c:2">(mol);
}

template void Rarey_smarts_part_210<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_210<RDKit::ROMol>(RDKit::ROMol &mol);
