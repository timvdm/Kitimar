#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_63(Mol &mol)
{
    // SMARTS 311 - 315
    validate<"[#6](-[#6]#[#7])(-[#6]#[#7])-[#6](-[$([#6]#[#7]),$([#6]=[#7])])-[#6]#[#7]">(mol);
    validate<"[#6](-[#6]#[#7])(-[#6]#[#7])=[#6](-[#16])-[#16]">(mol);
    validate<"[#6](-[#6]#[#7])(-[#6]#[#7])=[#6]-c:1:c:c:c:c:c:1">(mol);
    validate<"[#6](-[#6]#[#7])(-[#6]#[#7])=[#7]-[#7](-[#1])-c:1:c:c:c:c:c:1">(mol);
    validate<"[#6](-[#6]:[#6])(-[#6]:[#6])(-[#6]:[#6])-[#16]-[#6]:[#6]-[#6](=[#8])-[#8]-[#1]">(mol);
}

template void Rarey_smarts_part_63<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_63<RDKit::ROMol>(RDKit::ROMol &mol);
