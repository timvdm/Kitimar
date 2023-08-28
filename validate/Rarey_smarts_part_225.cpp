#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_225(Mol &mol)
{
    // SMARTS 1121 - 1125
    validate<"c:1:c:c(:c(:c:c:1)-[#6]=[#7]-[#7])-[#8]-[#1]">(mol);
    validate<"c:1:c:c(:c(:c:c:1)-[#8]-[#1])-[#8]-[#1]">(mol);
    validate<"c:1:c:c(:c:c-2:c:1-[#6](=[#6](-[#1])-[#6](=[#8])-[#8]-2)-c:3:c:c:c:c:c:3)-[#8]-[#6](-[#1])(-[#1])-[#6]:[#8]:[#6]">(mol);
    validate<"c:1:c:c(:c:c:c:1-[#6;X4]-c:2:c:c:c(:c:c:2)-[#7](-[$([#1]),$([#6;X4])])-[$([#1]),$([#6;X4])])-[#7](-[$([#1]),$([#6;X4])])-[$([#1]),$([#6;X4])]">(mol);
    validate<"c:1:c:c(:c:c:c:1-[#7](-[#1])-[#16](=[#8])=[#8])-[#7](-[#1])-[#16](=[#8])=[#8]">(mol);
}

template void Rarey_smarts_part_225<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_225<RDKit::ROMol>(RDKit::ROMol &mol);
