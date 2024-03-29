#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_68(Mol &mol)
{
    // SMARTS 336 - 340
    validate<"[#6]-1(=[!#6&!#1])-[#6](-[#7]=[#6]-[#16]-1)=[#8]">(mol);
    validate<"[#6]-1(=[#6](-!@[#6](=[#8])-[#7]-[#6](-[#1])-[#1])-[#16]-[#6](-[#7]-1-[$([#6](-[#1])(-[#1])-[#6](-[#1])=[#6](-[#1])-[#1]),$([#6]:[#6])])=[#16])-[$([#7]-[#6](=[#8])-[#6]:[#6]),$([#7](-[#1])-[#1])]">(mol);
    validate<"[#6]-1(=[#6](-!@[#6]=[#7])-[#16]-[#6](-[#7]-1)=[#8])-[$([F,Cl,Br,I]),$([#7+](:[#6]):[#6])]">(mol);
    validate<"[#6]-1(=[#6](-[#6](-[#1])(-[#6])-[#6])-[#16]-[#6](-[#7]-1-[$([#1]),$([#6](-[#1])-[#1])])=[#8])-[#16]-[#6;R]">(mol);
    validate<"[#6]-1(=[#6](-[#6](-[#6](-[#6](-[#6]-1(-[#1])-[#1])(-[#1])-[#6](=[#8])-[#6])(-[#1])-[#6](=[#8])-[#8]-[#1])(-[#1])-[#1])-[#6]:[#6])-[#6]:[#6]">(mol);
}

template void Rarey_smarts_part_68<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_68<RDKit::ROMol>(RDKit::ROMol &mol);
