#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_78(Mol &mol)
{
    // SMARTS 386 - 390
    validate<"[#6]:2(:[#6](-[#6](-[#1])-[#1]):[#6]-1:[#6](-[#7]=[#6](-[#7](-[#6]-1=[!#6&!#1;X1])-[#6](-[#1])-[$([#6](=[#8])-[#8]),$([#6]:[#6])])-[$([#1]),$([#16]-[#6](-[#1])-[#1])]):[!#6&!#1;X2]:2)-[#6](-[#1])(-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"[#6]:[#15$(a1aaaa1)]">(mol);
    validate<"[#6]:[#15$(a1aaaaa1)]">(mol);
    validate<"[#6]:[#16$(a1aaaa1)]">(mol);
    validate<"[#6]:[#16$(a1aaaaa1)]">(mol);
}

template void Rarey_smarts_part_78<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_78<RDKit::ROMol>(RDKit::ROMol &mol);
