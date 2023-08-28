#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_94(Mol &mol)
{
    // SMARTS 466 - 470
    validate<"[#7](-c:1:c:c:c:c:c:1)-c2[n+]c(cs2)-c:3:c:c:c:c:c:3">(mol);
    validate<"[#7]-1(-[#1])-[#6](=[#16])-[#6](-[#1])(-[#6]#[#7])-[#6](-[#1])(-[#6]:[#6])-[#6](=[#6]-1-[#6]:[#6])-[#1]">(mol);
    validate<"[#7]-1(-[#1])-[#7]=[#6](-[#7]-[#1])-[#16]-[#6](=[#6]-1-[#6]:[#6])-[#6]:[#6]">(mol);
    validate<"[#7]-1(-[#6](-[#1])-[#1])-[#6](=[#16])-[#7](-[#6]:[#6])-[#6](=[#7]-[#6]:[#6])-[#6]-1=[#7]-[#6]:[#6]">(mol);
    validate<"[#7]-1(-[$([#6;X4]),$([#1])])-[#6]=,:[#6](-[#6](=[#8])-[#6]:[#6]:[#6])-[#6](-[#6])-[#6](=[#6]-1-[#6](-[#1])(-[#1])-[#1])-[$([#6]=[#8]),$([#6]#[#7])]">(mol);
}

template void Rarey_smarts_part_94<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_94<RDKit::ROMol>(RDKit::ROMol &mol);
