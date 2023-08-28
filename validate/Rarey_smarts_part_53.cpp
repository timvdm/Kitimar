#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_53(Mol &mol)
{
    // SMARTS 261 - 265
    validate<"[#6;X4]-1-[#6](=[#8])-[#7]-[#7]-[#6]-1=[#8]">(mol);
    validate<"[#6;X4]-[#16;X2]-[#6](=[#7]-[!#1]:[!#1]:[!#1]:[!#1])-[#7](-[#1])-[#7]=[#6]">(mol);
    validate<"[#6;X4]-[#7+](-[#6;X4]-[#8]-[#1])=[#6]-[#16]-[#6]-[#1]">(mol);
    validate<"[#6;X4]-[#7](-[#1])-[#6](-[#6]:[#6])=[#6](-[#1])-[#6](=[#16])-[#7](-[#1])-c:1:c:c:c:c:c:1">(mol);
    validate<"[#6;X4]-[#7](-[#6;X4])-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#6]2=,:[#7][#6]:[#6]:[!#1]2)-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_53<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_53<RDKit::ROMol>(RDKit::ROMol &mol);
