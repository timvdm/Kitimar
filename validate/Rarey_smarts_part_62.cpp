#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_62(Mol &mol)
{
    // SMARTS 306 - 310
    validate<"[#6](-[#1])-[#7](-[#1])-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1])-[#7](-[#1])-[#6]-[#1]">(mol);
    validate<"[#6](-[#1])-[#7](-[#1])-c:1:c(:c(:c(:s:1)-[#6]-[#1])-[#6]-[#1])-[#6](=[#8])-[#7](-[#1])-[#6]:[#6]">(mol);
    validate<"[#6](-[#1])-[#7](-[#1])-c:1:n:c(:c:s:1)-c2cnc3n2ccs3">(mol);
    validate<"[#6](-[#1])-[#7](-[#6](-[#1])-[#1])-c:1:c(:c(:c(:c(:c:1-[#1])-[$([#1]),$([#6](-[#1])-[#1])])-[#6](-[#1])-[$([#1]),$([#6]-[#1])])-[#1])-[#1]">(mol);
    validate<"[#6](-[#6]#[#7])(-[#6]#[#7])-[#6](-[#7](-[#1])-[#1])=[#6]-[#6]#[#7]">(mol);
}

template void Rarey_smarts_part_62<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_62<RDKit::ROMol>(RDKit::ROMol &mol);
