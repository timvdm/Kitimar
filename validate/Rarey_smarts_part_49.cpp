#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_49(Mol &mol)
{
    // SMARTS 241 - 245
    validate<"[#16]-1-[#6](=[#7]-[#6]:[#6])-[#7](-[$([#1]),$([#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#8]),$([#6]:[#6])])-[#6](=[#8])-[#6]-1=[#6](-[#1])-[$([#6]:[#6]:[#6]-[#17]),$([#6]:[!#6&!#1])]">(mol);
    validate<"[#16]-1-[#6](=[#7]-[#7]-[#1])-[#16]-[#6](=[#7]-[#6]:[#6])-[#6]-1=[#7]-[#6]:[#6]">(mol);
    validate<"[#16]-1-[#6](=[#8])-[#7]-[#6](=[#16])-[#6]-1=[#6](-[#1])-[#6]:[#6]">(mol);
    validate<"[#16]-1-[#6](=[#8])-[#7]-[#6](=[#8])-[#6]-1=[#6](-[#1])-[$([#6]-[#35]),$([#6]:[#6](-[#1]):[#6](-[F,Cl,Br,I]):[#6]:[#6]-[F,Cl,Br,I]),$([#6]:[#6](-[#1]):[#6](-[#1]):[#6]-[#16]-[#6](-[#1])-[#1]),$([#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#8]-[#6](-[#1])-[#1]),$([#6]:1:[#6](-[#6](-[#1])-[#1]):[#7](-[#6](-[#1])-[#1]):[#6](-[#6](-[#1])-[#1]):[#6]:1)]">(mol);
    validate<"[#16]:[#15$(a1aaaa1)]">(mol);
}

template void Rarey_smarts_part_49<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_49<RDKit::ROMol>(RDKit::ROMol &mol);
