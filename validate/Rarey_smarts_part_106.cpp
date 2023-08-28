#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_106(Mol &mol)
{
    // SMARTS 526 - 530
    validate<"[#8]=[#6]-!@n:1:c:c:c-2:c:1-[#7](-[#1])-[#6](=[#16])-[#7]-2-[#1]">(mol);
    validate<"[#8]=[#6]-1-[#6;X4]-[#6]-[#6](=[#8])-c:2:c:c:c:c:c-1:2">(mol);
    validate<"[#8]=[#6]-1-[#6](=[#6]-[#6](=[#7]-[#7]-1)-[#6]=[#8])-[#6]#[#7]">(mol);
    validate<"[#8]=[#6]-1-[#6](=[#7]-[#7]-[#6]-[#6]-1)-[#6]#[#7]">(mol);
    validate<"[#8]=[#6]-1-[#6]:[#6]-[#6](-[#1])(-[#1])-[#7]-[#6]-1=[#6]-[#1]">(mol);
}

template void Rarey_smarts_part_106<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_106<RDKit::ROMol>(RDKit::ROMol &mol);
