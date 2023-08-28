#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_232(Mol &mol)
{
    // SMARTS 1156 - 1160
    validate<"c:2(:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#6]=[#6]-[#6]-3=[#7])-[#7]">(mol);
    validate<"c:2(:c:1:c:c:c:c:c:1:c-3:c(:c:2)-[#6](-c:4:c:c:c:c:c-3:4)=[#8])-[#8]-[#1]">(mol);
    validate<"c:2(:c:1:c:c:c:c:c:1:n:n:c:2)-[#6](-[#6]:[#6])-[#6]#[#7]">(mol);
    validate<"c:2(:n:c:1:c(:c(:c:c(:c:1-[#1])-[F,Cl,Br,I])-[#1]):n:2-[#1])-[#16]-[#6](-[#1])(-[#1])-[#6](=[#8])-[#7](-[#1])-[#6]:[#6]">(mol);
    validate<"c:2-3:c(:c:c:1:c:c:c:c:c:1:c:2)-[#7](-[#6](-[#1])-[#1])-[#6](=[#8])-[#6](=[#7]-3)-[#6]:[#6]-[#7](-[#1])-[#6](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_232<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_232<RDKit::ROMol>(RDKit::ROMol &mol);
