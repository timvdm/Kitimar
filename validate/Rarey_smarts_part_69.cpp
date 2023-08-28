#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_69(Mol &mol)
{
    // SMARTS 341 - 345
    validate<"[#6]-1(=[#6](-[#6](=[#8])-[#7]-[#6](=[#7]-1)-[!#6&!#1])-[#6]#[#7])-[#6]">(mol);
    validate<"[#6]-1(=[#6](-[$([#1]),$([#6](-[#1])-[#1]),$([#6]=[#8])])-[#16]-[#6](-[#7]-1-[$([#1]),$([#6]-[#1]),$([#6]:[#6])])=[#7;!R])-[$([#6](-[#1])-[#1]),$([#6]:[#6])]">(mol);
    validate<"[#6]-1(=[#6])-[#6](-[#7]=[#6]-[#16]-1)=[#8]">(mol);
    validate<"[#6]-1(=[#6])-[#6](=[#8])-[#7]-[#7]-[#6]-1=[#8]">(mol);
    validate<"[#6]-1(=[#6])-[#6]=[#7]-[!#6&!#1]-[#6]-1=[#8]">(mol);
}

template void Rarey_smarts_part_69<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_69<RDKit::ROMol>(RDKit::ROMol &mol);
