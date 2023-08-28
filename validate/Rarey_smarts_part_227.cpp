#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_227(Mol &mol)
{
    // SMARTS 1131 - 1135
    validate<"c:1:c:c-2:c(:c:c:1)-[#16]-c3c(-[#7]-2)cc(s3)-[#6](-[#1])-[#1]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#6](-c3cccc4noc-2c34)=[#8]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#6](-c:3:c(-[$([#16;X2]),$([#6;X4])]-2):c:c:c(:c:3)-[$([#1]),$([#17]),$([#6;X4])])=[#6]-[#6]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#6](=[#6](-[#6]-2=[#8])-[#6])-[#8]-[#1]">(mol);
    validate<"c:1:c:c-2:c(:c:c:1)-[#6]-3-[#6](-[#6]-[#7]-2)-[#6]-[#6]=[#6]-3">(mol);
}

template void Rarey_smarts_part_227<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_227<RDKit::ROMol>(RDKit::ROMol &mol);
