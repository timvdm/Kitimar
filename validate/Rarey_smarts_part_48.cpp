#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_48(Mol &mol)
{
    // SMARTS 236 - 240
    validate<"[#16X2][OX2H,OX1H0-]">(mol);
    validate<"[#16X2][OX2H0]">(mol);
    validate<"[#16](=[#8])(=[#8])(-c:1:c:n(-[#6](-[#1])-[#1]):c:n:1)-[#7](-[#1])-c:2:c:n(:n:c:2)-[#6](-[#1])(-[#1])-[#6]:[#6]-[#8]-[#6](-[#1])-[#1]">(mol);
    validate<"[#16](=[#8])(=[#8])-[#7](-[#1])-c:1:c(:c(:c(:s:1)-[#6]-[#1])-[#6]-[#1])-[#6](=[#8])-[#7]-[#1]">(mol);
    validate<"[#16]-1-[#6](=!@[#7]-[$([#1]),$([#7](-[#1])-[#6]:[#6])])-[#7](-[$([#1]),$([#6]:[#7]:[#6]:[#6]:[#16])])-[#6](=[#8])-[#6]-1=[#6](-[#1])-[#6]:[#6]-[$([#17]),$([#8]-[#6]-[#1])]">(mol);
}

template void Rarey_smarts_part_48<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_48<RDKit::ROMol>(RDKit::ROMol &mol);
