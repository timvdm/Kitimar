#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_211(Mol &mol)
{
    // SMARTS 1051 - 1055
    validate<"c:1(:c:c:c:c:c:1)-[#7](-[#1])-[#6](=[#16])-[#7]-[#7](-[#1])-[#6](-[#1])=[#6](-[#1])-[#6]=[#8]">(mol);
    validate<"c:1(:c:c:c:c:c:1)-[#7](-[#1])-[#6](=[#16])-[#7]-[#7](-[#1])-[#6](=[#8])-[#6]-2:[!#1]:[!#6&!#1]:[#6]:[#6]-2">(mol);
    validate<"c:1(:c:c:c:c:c:1)-[#7](-[#1])-[#6](=[#16])-[#7]-[#7](-[#1])-[#6]([#7;R])[#7;R]">(mol);
    validate<"c:1(:c:c:c:c:c:1)-[#7](-[#1])-[#6](=[#16])-[#7]-[#7](-[#1])-c:2:c:c:c:c:c:2">(mol);
    validate<"c:1(:c:c:c:c:c:1-[#7](-[#1])-[#7]=[#6])-[#6](=[#8])-[#8]-[#1]">(mol);
}

template void Rarey_smarts_part_211<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_211<RDKit::ROMol>(RDKit::ROMol &mol);
