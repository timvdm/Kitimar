#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_79(Mol &mol)
{
    // SMARTS 391 - 395
    validate<"[#6]:[#6$(a1aaaa1)]">(mol);
    validate<"[#6]:[#6$(a1aaaaa1)]">(mol);
    validate<"[#6]:[#6]-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#6]=[#8])-[#7]-2-[#6](=[#8])-[#6]-1(-[#1])-[#6](-[#1])(-[#1])-[#6]=[#6]-[#6](-[#1])(-[#1])-[#6]-1(-[#1])-[#6]-2=[#8]">(mol);
    validate<"[#6]:[#6]-[#6](-[#1])=[#6](-[#1])-[#6](-[#1])=[#7]-[#7](-[#6;X4])-[#6;X4]">(mol);
    validate<"[#6]:[#6]-[#6](-[#1])=[#6](-[#1])-[#6](-[#1])=[#7]-[#7]=[#6]">(mol);
}

template void Rarey_smarts_part_79<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_79<RDKit::ROMol>(RDKit::ROMol &mol);
