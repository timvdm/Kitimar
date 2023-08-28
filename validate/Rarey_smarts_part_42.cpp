#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_42(Mol &mol)
{
    // SMARTS 206 - 210
    validate<"[!#1]:[#6]-[#6](=[#16])-[#7](-[#1])-[#7](-[#1])-[#6]:[!#1]">(mol);
    validate<"[!#1]:[#6]-[#6]-1=[#6](-[#1])-[#6](=[#6](-[#6]#[#7])-[#6](=[#8])-[#7]-1-[#1])-[#6]:[#8]">(mol);
    validate<"[!#6&!#1]=[#6]-1-[#6]=,:[#6]-[#6](=[!#6&!#1])-[#6]=,:[#6]-1">(mol);
    validate<"[!#6]-[CH3]">(mol);
    validate<"[!#6]~*~*~*~[D3]">(mol);
}

template void Rarey_smarts_part_42<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_42<RDKit::ROMol>(RDKit::ROMol &mol);
