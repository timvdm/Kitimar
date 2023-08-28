#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_241(Mol &mol)
{
    // SMARTS 1201 - 1205
    validate<"n2(-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1])-[#6](=[#8])-[#7](-[#1])-[#6](-[#1])(-[#6](-[#1])-[#1])-[#6](-[#1])(-[#1])-[#8]-[#6]:[#6])c(c(-[#1])c(c2-[#1])-[#1])-[#1]">(mol);
    validate<"n2(-c:1:c(:c:c(:c(:c:1)-[#1])-[$([#7](-[#1])-[#1]),$([#6]:[#7])])-[#1])c(c(-[#1])c(c2-[#1])-[#1])-[#1]">(mol);
    validate<"n2(-c:1:c:c:c:c:c:1)c(c(-[#1])c(c2-[#6]=[#7]-[#8]-[#1])-[#1])-[#1]">(mol);
    validate<"n2c1ccccn1c(c2-[$([#6](-[!#1])=[#6](-[#1])-[#6]:[#6]),$([#6]:[#8]:[#6])])-[#7]-[#6]:[#6]">(mol);
    validate<"n2nc(c1cccc1c2-[#6])-[#6]">(mol);
}

template void Rarey_smarts_part_241<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_241<RDKit::ROMol>(RDKit::ROMol &mol);
