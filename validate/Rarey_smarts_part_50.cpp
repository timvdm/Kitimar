#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_50(Mol &mol)
{
    // SMARTS 246 - 250
    validate<"[#16]:[#15$(a1aaaaa1)]">(mol);
    validate<"[#16]:[#16$(a1aaaa1)]">(mol);
    validate<"[#16]:[#16$(a1aaaaa1)]">(mol);
    validate<"[#16]=[#6]-1-[#6]=,:[#6]-[!#6&!#1]-[#6]=,:[#6]-1">(mol);
    validate<"[#16]=[#6]-1-[#7](-[#1])-[#6]=[#6]-[#6]-2=[#6]-1-[#6](=[#8])-[#8]-[#6]-2=[#6]-[#1]">(mol);
}

template void Rarey_smarts_part_50<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_50<RDKit::ROMol>(RDKit::ROMol &mol);
