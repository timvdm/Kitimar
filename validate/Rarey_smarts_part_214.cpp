#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_214(Mol &mol)
{
    // SMARTS 1066 - 1070
    validate<"c:1-2:c(:c(:c(:c(:c:1-[#8]-[#6](-[#1])(-[#1])-[#8]-2)-[#6](-[#1])(-[#1])-[#7]-3-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#6]:[#6]-3)-[#1])-[#1])-[#1]">(mol);
    validate<"c:1-2:c(:c:c(:c:c:1-[F,Cl,Br,I])-[F,Cl,Br,I])-[#6](=[#6](-[#6](=[#8])-[#7](-[#1])-[#1])-[#6](=[#7]-[#1])-[#8]-2)-[#1]">(mol);
    validate<"c:1-2:c(:c:c:c(:c:1-[#8]-[#6](-[#1])(-[#1])-[#7](-[#6]:[#6]-[#8]-[#6](-[#1])-[#1])-[#6]-2(-[#1])-[#1])-[#1])-[#1]">(mol);
    validate<"c:1-2:c(:c:c:c:c:1)-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#7]=[#6]-2-[#16;X2]-[#6](-[#1])(-[#1])-[#6](=[#8])-c:3:c:c:c:c:c:3">(mol);
    validate<"c:1-2:c(:c:c:c:c:1)-[#6](=[#8])-[#6;X4]-[#6]-2=[#8]">(mol);
}

template void Rarey_smarts_part_214<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_214<RDKit::ROMol>(RDKit::ROMol &mol);
