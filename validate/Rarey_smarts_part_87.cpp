#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_87(Mol &mol)
{
    // SMARTS 431 - 435
    validate<"[#7;!R]=[#7]">(mol);
    validate<"[#7;X2v4+0]">(mol);
    validate<"[#7]">(mol);
    validate<"[#7](-[#1])(-[#1])-[#6]-1=[#6](-[#6]#[#7])-[#6](-[#1])(-[#6]:[#6])-[#16]-[#6;X4]-[#16]-1">(mol);
    validate<"[#7](-[#1])(-[#1])-[#6]-1=[#6](-[#6]#[#7])-[#6](-[#1])(-[#6]:[#6])-[#6](=[#6](-[#6]:[#6])-[#8]-1)-[#6]#[#7]">(mol);
}

template void Rarey_smarts_part_87<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_87<RDKit::ROMol>(RDKit::ROMol &mol);
