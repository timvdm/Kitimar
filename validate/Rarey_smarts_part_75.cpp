#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_75(Mol &mol)
{
    // SMARTS 371 - 375
    validate<"[#6]-3(=[#8])-[#6](=[#6](-[#1])-[#7](-[#1])-c:1:c:c:c:c:c:1-[#6](=[#8])-[#8]-[#1])-[#7]=[#6](-c:2:c:c:c:c:c:2)-[#8]-3">(mol);
    validate<"[#6]-[#16;X2]-c:1:n:c(:c:s:1)-[#1]">(mol);
    validate<"[#6]-[#6](=[!#6&!#1;!R])-[#6](=[!#6&!#1;!R])-[$([#6]),$([#16](=[#8])=[#8])]">(mol);
    validate<"[#6]-[#6](=[#16])-[#1]">(mol);
    validate<"[#6]-[#6](=[#16])-[#6]">(mol);
}

template void Rarey_smarts_part_75<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_75<RDKit::ROMol>(RDKit::ROMol &mol);
