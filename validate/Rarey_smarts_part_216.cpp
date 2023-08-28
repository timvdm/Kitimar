#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_216(Mol &mol)
{
    // SMARTS 1076 - 1080
    validate<"c:1-3:c(:c:2:c(:c:c:1-[Br]):o:c:c:2)-[#6](=[#6]-[#6](=[#8])-[#8]-3)-[#1]">(mol);
    validate<"c:1-3:c(:c:c:c:c:1)-[#16]-[#6](=[#7]-[#7]=[#6]-2-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-2)-[#7]-3-[#6](-[#1])-[#1]">(mol);
    validate<"c:1-3:c(:c:c:c:c:1)-[#6](=[#6](-[#6](=[#8])-[#7](-[#1])-c:2:n:c(:c:s:2)-[#6]:[#16]:[#6]-[#1])-[#6](=[#8])-[#8]-3)-[#1]">(mol);
    validate<"c:1-3:c(:c:c:c:c:1)-[#6](=[#6](-[#6](=[#8])-[#7](-[#1])-c:2:n:o:c:c:2-[Br])-[#6](=[#8])-[#8]-3)-[#1]">(mol);
    validate<"c:1-3:c:2:c(:c(:c:c:1)-[#7]):c:c:c:c:2-[#6](-[#6]=[#6]-3-[#6](-[F])(-[F])-[F])=[#8]">(mol);
}

template void Rarey_smarts_part_216<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_216<RDKit::ROMol>(RDKit::ROMol &mol);
