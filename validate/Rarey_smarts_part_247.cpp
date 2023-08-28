#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_247(Mol &mol)
{
    // SMARTS 1231 - 1235
    validate<"s2c:1:n:c:n:c(:c:1c(c2-[#6](-[#1])-[#1])-[#6](-[#1])-[#1])-[#7]-[#7]=[#6]-c3ccco3">(mol);
    validate<"s:1:c(:[n+](-[#6](-[#1])-[#1]):c(:c:1-[#1])-[#6])-[#7](-[#1])-c:2:c:c:c:c:c:2[$([#6](-[#1])-[#1]),$([#6]:[#6])]">(mol);
    validate<"s:1:c(:c(-[#1]):c(:c:1-[#6](=[#8])-[#7](-[#1])-[#7]-[#1])-[#8]-[#6](-[#1])-[#1])-[#1]">(mol);
    validate<"s:1:c(:c(-[#1]):c(:c:1-[#6]-3=[#7]-c:2:c:c:c:c:c:2-[#6](=[#7]-[#7]-3-[#1])-c:4:c:c:n:c:c:4)-[#1])-[#1]">(mol);
    validate<"s:1:c:c:c(:c:1-[#1])-c:2:c:s:c(:n:2)-[#7](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_247<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_247<RDKit::ROMol>(RDKit::ROMol &mol);
