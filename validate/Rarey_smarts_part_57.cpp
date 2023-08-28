#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_57(Mol &mol)
{
    // SMARTS 281 - 285
    validate<"[#6](-[#1])(-[#1])-[#6](-[#1])(-[#6]#[#7])-[#6](=[#8])-[#6]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#6](-[#8]-[#1])=[#6](-[#6](=[#8])-[#6](-[#1])-[#1])-[#6](-[#1])-[#6]#[#6]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#1])-[#6]=[#7]-[#7](-[#1])-c1nc(c(-[#1])s1)-[#6]:[#6]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-[#6](-[#1])=[#6]-[#6](=[#8])-c:1:c(-[#16;X2]):s:c(:c:1)-[$([#6]#[#7]),$([#6]=[#8])]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])-[#1])-[#6]-2=[#6](-[#1])-c:1:c(:c:c:c:c:1)-[#16;X2]-c:3:c-2:c:c:c:c:3">(mol);
}

template void Rarey_smarts_part_57<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_57<RDKit::ROMol>(RDKit::ROMol &mol);
