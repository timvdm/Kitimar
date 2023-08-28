#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_245(Mol &mol)
{
    // SMARTS 1221 - 1225
    validate<"n@a[#8]">(mol);
    validate<"o:1:c(:c(-[#1]):c(:c:1-[#6](-[#1])(-[#1])-[#7](-[#1])-[#6](=[#16])-[#7](-[#6]-[#1])-[#6](-[#1])(-[#1])-c:2:c:c:c:c:c:2)-[#1])-[#1]">(mol);
    validate<"o:1:c(:c:c:2:c:1:c(:c(:c(:c:2-[#1])-[#8]-[#6](-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])-[#1])-[#6](~[#8])~[#8]">(mol);
    validate<"p-[Br]">(mol);
    validate<"p-[Cl]">(mol);
}

template void Rarey_smarts_part_245<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_245<RDKit::ROMol>(RDKit::ROMol &mol);
