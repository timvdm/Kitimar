#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_240(Mol &mol)
{
    // SMARTS 1196 - 1200
    validate<"n2(-[#6](-[#1])-c:1:c(:c(:c:c(:c:1-[#1])-[#1])-[#1])-[#1])c(c(-[#1])c(c2-[#6]-[#1])-[#1])-[#6]-[#1]">(mol);
    validate<"n2(-[#6]:1:[!#1]:[!#6&!#1]:[!#1]:[#6]:1-[#1])c(c(-[#1])c(c2-[#6;X4])-[#1])-[#6;X4]">(mol);
    validate<"n2(-[#6]:1:[!#1]:[#6]:[#6]:[#6]:[#6]:1)c(cc(c2-[#6;X4])-[#1])-[#6;X4]">(mol);
    validate<"n2(-[#6]:1:[!#1]:[#6]:[#6]:[#6]:[#6]:1)c(cc(c2-[#6]:[#6])-[#1])-[#6;X4]">(mol);
    validate<"n2(-[#6]:1:[#6](-[#6]#[#7]):[#6]:[#6]:[!#6&!#1]:1)c(c(-[#1])c(c2)-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_240<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_240<RDKit::ROMol>(RDKit::ROMol &mol);
