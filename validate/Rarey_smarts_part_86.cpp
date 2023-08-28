#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_86(Mol &mol)
{
    // SMARTS 426 - 430
    validate<"[#7+](:[!#1]:[!#1]:[!#1])-[!#1]=[#8]">(mol);
    validate<"[#7+]([#6]:[#6])=,:[#6]-[#6](-[#1])=[#6]-[#7](-[#6;X4])-[#6]">(mol);
    validate<"[#7+]:1(:[#6]:[#6]:[!#1]:c:2:c:1:c(:c(-[$([#1]),$([#7])]):c:c:2)-[#1])-[$([#6](-[#1])(-[#1])-[#1]),$([#8;X1]),$([#6](-[#1])(-[#1])-[#6](-[#1])=[#6](-[#1])-[#1]),$([#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#8]-[#1]),$([#6](-[#1])(-[#1])-[#6](=[#8])-[#6]),$([#6](-[#1])(-[#1])-[#6](=[#8])-[#7](-[#1])-[#6]:[#6]),$([#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#1])]">(mol);
    validate<"[#7,#8,#16;-0,-1;!$([o,s,nX3]);!$([Nv5,Pv4,Pv5,Sv4,Sv6])]">(mol);
    validate<"[#7;!R]=[#6]-2-[#6](=[#8])-c:1:c:c:c:c:c:1-[#16]-2">(mol);
}

template void Rarey_smarts_part_86<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_86<RDKit::ROMol>(RDKit::ROMol &mol);
