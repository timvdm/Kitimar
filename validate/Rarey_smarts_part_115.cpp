#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_115(Mol &mol)
{
    // SMARTS 571 - 575
    validate<"[$([#1]),$([#6](-[#1])-[#1]),$([#6]:[#6])]-c:1:c(:c(:c(:s:1)-[#7](-[#1])-[#6](=[#8])-[#6])-[#6](=[#8])-[#8])-[$([#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1),$([#6]:1:[#16]:[#6]:[#6]:[#6]:1)]">(mol);
    validate<"[$([#1]),$([#6](-[#1])-[#1])]-[#8]-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1])-[#7](-[#1])-[#6](-[#1])(-[#1])-c:2:n:c:c:n:2">(mol);
    validate<"[$([#6+0]);!$(C(F)(F)F);!$(c(:[!c]):[!c])!$([#6]=,#[!#6])]">(mol);
    //validate<"[$([#6X4@](*)(*)(*)*),$([#6X4@H](*)(*)*)]">(mol); // FIXME: stereo
    validate<"[$([#6]-!@[A!C!H]-!@[#6])]">(mol);
}

template void Rarey_smarts_part_115<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_115<RDKit::ROMol>(RDKit::ROMol &mol);
