#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_67(Mol &mol)
{
    // SMARTS 331 - 335
    validate<"[#6]-1(-[#6](=[#8])-[#7]-[#6](=[#8])-[#7]-[#6]-1=[#8])=[#7]">(mol);
    validate<"[#6]-1(-[#6](~[!#6&!#1]~[#6]-[!#6&!#1]-[#6]-1=[!#6&!#1])~[!#6&!#1])=[#6;!R]-[#1]">(mol);
    validate<"[#6]-1(-[#6]=,:[#6]-[#6]=,:[#6]-[#6]-1=[!#6&!#1])=[!#6&!#1]">(mol);
    validate<"[#6]-1(-[#6]=[#8])(-[#6]:[#6])-[#16;X2]-[#6]=[#7]-[#7]-1-[#1]">(mol);
    validate<"[#6]-1(=[!#1]-[!#1]=[!#1]-[#7](-[#6]-1=[#16])-[#1])-[#6]#[#7]">(mol);
}

template void Rarey_smarts_part_67<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_67<RDKit::ROMol>(RDKit::ROMol &mol);
