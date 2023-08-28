#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_100(Mol &mol)
{
    // SMARTS 496 - 500
    validate<"[#7]-4(-c:1:c:c:c:c:c:1)-[#6](=[#8])-[#16]-[#6](-[#1])(-[#7](-[#1])-c:2:c:c:c:c:3:c:c:c:c:c:2:3)-[#6]-4=[#8]">(mol);
    validate<"[#7]-[#6;X4]-c:1:c:c:c:c:c:1-[#8]-[#1]">(mol);
    validate<"[#7]-[#6]=!@[#6]-2-[#6](=[#8])-c:1:c:c:c:c:c:1-[!#6&!#1]-2">(mol);
    validate<"[#7]:[#15$(a1aaaa1)]">(mol);
    validate<"[#7]:[#16$(a1aaaaa1)]">(mol);
}

template void Rarey_smarts_part_100<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_100<RDKit::ROMol>(RDKit::ROMol &mol);
