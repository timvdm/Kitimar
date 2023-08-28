#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_238(Mol &mol)
{
    // SMARTS 1186 - 1190
    validate<"n1(-[#6;X4])c(c(-[#1])c(c1-[#6]:[#6])-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"n1(-[#6](-[#1])-[#1])c(c(-[#6](=[#8])-[#6])c(c1-[#6]:[#6])-[#6])-[#6](-[#1])-[#1]">(mol);
    validate<"n1(-[#6])c(c(-[#1])c(c1-[#6](-[#1])(-[#1])-[#7](-[#1])-[#6](=[#16])-[#7]-[#1])-[#1])-[#1]">(mol);
    validate<"n1(-[#6])c(c(-[#1])c(c1-[#6](-[#1])=[#6](-[#6]#[#7])-c:2:n:c:c:s:2)-[#1])-[#1]">(mol);
    validate<"n1(-[#6])c(c(-[#1])c(c1-[#6](-[#1])=[#6]-2-[#6](=[#8])-[!#6&!#1]-[#6]=,:[!#1]-2)-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_238<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_238<RDKit::ROMol>(RDKit::ROMol &mol);
