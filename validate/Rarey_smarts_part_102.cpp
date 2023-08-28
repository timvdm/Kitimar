#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_102(Mol &mol)
{
    // SMARTS 506 - 510
    validate<"[#7]=[#6]-1-[#7](-[#1])-[#6](=[#6](-[#7]-[#1])-[#7]=[#7]-1)-[#7]-[#1]">(mol);
    validate<"[#7]=[#6]-1-[#7]=[#6]-[#7]-[#16]-1">(mol);
    validate<"[#7]~*(~*)~*">(mol);
    validate<"[#7]~*~*~*~*~*~*~*~*~[#8]">(mol);
    validate<"[#7]~[#6]:1:[#7]:[#7]:[#6](:[$([#7]),$([#6]-[#1]),$([#6]-[#7]-[#1])]:[$([#7]),$([#6]-[#7])]:1)-[$([#7]-[#1]),$([#8]-[#6](-[#1])-[#1])]">(mol);
}

template void Rarey_smarts_part_102<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_102<RDKit::ROMol>(RDKit::ROMol &mol);
