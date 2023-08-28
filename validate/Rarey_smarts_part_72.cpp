#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_72(Mol &mol)
{
    // SMARTS 356 - 360
    validate<"[#6]-1=[!#1]-[!#6&!#1]-[#6](-[#6]-1=[!#6&!#1;!R])=[#8]">(mol);
    validate<"[#6]-1=[#6](-[#16]-[#6](-[#6]=[#6]-1)=[#16])-[#7]">(mol);
    validate<"[#6]-1=[#6]-[#6](-[#8]-[#6]-1-[#8])(-[#8])-[#6]">(mol);
    validate<"[#6]-1=[#6]-[#7](-[#6](-c:2:c-1:c:c:c:c:2)(-[#6]#[#7])-[#6](=[#16])-[#16])-[#6]=[#8]">(mol);
    validate<"[#6]-1=[#6]-[#7]-[#6](-[#16]-[#6;X4]-1)=[#16]">(mol);
}

template void Rarey_smarts_part_72<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_72<RDKit::ROMol>(RDKit::ROMol &mol);
