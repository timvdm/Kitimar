#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_95(Mol &mol)
{
    // SMARTS 471 - 475
    validate<"[#7]-1-2-[#6](=[#7]-[#6](=[#8])-[#6](=[#7]-1)-[#6](-[#1])-[#1])-[#16]-[#6](=[#6](-[#1])-[#6]:[#6])-[#6]-2=[#8]">(mol);
    validate<"[#7]-1-[#6](=[#16])-[#16]-[#6;X4]-[#6]-1=[#8]">(mol);
    validate<"[#7]-1-[#6](=[#16])-[#16]-[#6](=[#6])-[#6]-1=[#8]">(mol);
    validate<"[#7]-1-[#6](=[#8])-[#6](=[#6](-[#6])-[#16]-[#6]-1=[#16])-[#1]">(mol);
    validate<"[#7]-1=[#6](-[#7](-[#6](-[#6](-[#6]-1(-[#1])-[#6]:[#6])(-[#1])-[#1])=[#8])-[#1])-[#7]-[#1]">(mol);
}

template void Rarey_smarts_part_95<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_95<RDKit::ROMol>(RDKit::ROMol &mol);
