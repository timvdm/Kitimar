#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_65(Mol &mol)
{
    // SMARTS 321 - 325
    validate<"[#6](=[#6](-[#1])-[#6](-[#1])(-[#1])-n:1:c(:n:c(:c:1-[#1])-[#1])-[#1])(-[#6]:[#6])-[#6]:[#6]">(mol);
    validate<"[#6](=[#8])(-[#7]-1-[#6]-[#6]-[#16]-[#6]-[#6]-1)-c:2:c(:c(:c(:c(:c:2-[#16]-[#6](-[#1])-[#1])-[#1])-[#1])-[#1])-[#1]">(mol);
    validate<"[#6](=[#8])-[#6](-[#1])=[#6](-[#8]-[#1])-[#6](-[#8]-[#1])=[#6](-[#1])-[#6](=[#8])-[#6]">(mol);
    validate<"[#6](=[#8])-[#6]-1=[#6]-[#7]-c:2:c(-[#16]-1):c:c:c:c:2">(mol);
    validate<"[#6]-1(-[#6]#[#7])(-[#6]#[#7])-[#6](-[#1])(-[#6](=[#8])-[#6])-[#6]-1-[#1]">(mol);
}

template void Rarey_smarts_part_65<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_65<RDKit::ROMol>(RDKit::ROMol &mol);
