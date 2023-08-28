#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_41(Mol &mol)
{
    // SMARTS 201 - 205
    validate<"[!#1]">(mol);
    validate<"[!#1]:1:[!#1]-2:[!#1](:[!#1]:[!#1]:[!#1]:1)-[#7](-[#1])-[#7](-[#6]-2=[#8])-[#6]">(mol);
    validate<"[!#1]:1:[!#1]:[!#1]:[!#1](:[!#1]:[!#1]:1)-[#6](-[#1])=[#6](-[#1])-[#6](-[#7](-[#1])-[#7](-[#1])-c2nnnn2-[#6])=[#8]">(mol);
    validate<"[!#1]:1:[!#1]:[!#1]:[!#1](:[!#1]:[!#1]:1)-[#6](-[#1])=[#6](-[#1])-[#6](-[#7]-c:2:c:c:c:3:c(:c:2):c:c:c(:n:3)-[#7](-[#6])-[#6])=[#8]">(mol);
    validate<"[!#1]:[!#1]-[#6](-[$([#1]),$([#6]#[#7])])=[#6]-1-[#6]=,:[#6]-[#6](=[$([#8]),$([#7;!R])])-[#6]=,:[#6]-1">(mol);
}

template void Rarey_smarts_part_41<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_41<RDKit::ROMol>(RDKit::ROMol &mol);
