#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_239(Mol &mol)
{
    // SMARTS 1191 - 1195
    validate<"n1(-[#6])c(c(-[#1])c(c1-[#6]=[#7]-[#7])-[#1])-[#1]">(mol);
    validate<"n1-2cccc1-[#6]=[#7](-[#6])-[#6]-[#6]-2">(mol);
    validate<"n1nnnc2cccc12">(mol);
    validate<"n1nscc1-c2nc(no2)-[#6]:[#6]">(mol);
    validate<"n2(-[#6](-[#1])-[#1])c-1c(-[#6]:[#6]-[#6]-1=[#8])cc2-[#6](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_239<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_239<RDKit::ROMol>(RDKit::ROMol &mol);
