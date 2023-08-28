#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_244(Mol &mol)
{
    // SMARTS 1216 - 1220
    validate<"n:1:c(:n(:c(:c:1-c:2:c:c:c:c:c:2)-c:3:c:c:c:c:c:3)-[#7]=!@[#6])-[#7](-[#1])-[#1]">(mol);
    validate<"n:1:c(:n(:c:2:c:1:c:c:c:c:2)-[#6](-[#1])-[#1])-[#16]-[#6](-[#1])(-[#1])-[#6](=[#8])-[#7](-[#1])-[#7]=[#6](-[#1])-[#6](-[#1])=[#6]-[#1]">(mol);
    validate<"n:1:c(:n:c(:n:c:1-[#7](-[#6](-[#1])-[#1])-[#6](-[#1])-[#1])-[#7](-[#6](-[#1])-[#1])-[#6](-[#1])-[#1])-[#7](-[#6]-[#1])-[#6]=[#8]">(mol);
    validate<"n:1:c3:c(:c:c2:c:1nc(s2)-[#7])sc(n3)-[#7]">(mol);
    validate<"n:1:c:c:c(:c:1-[#6](-[#1])-[#1])-[#6](-[#1])=[#6]-2-[#6](=[#8])-[#7]-[#6](=[!#6&!#1])-[#7]-2">(mol);
}

template void Rarey_smarts_part_244<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_244<RDKit::ROMol>(RDKit::ROMol &mol);
