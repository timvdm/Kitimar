#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_70(Mol &mol)
{
    // SMARTS 346 - 350
    validate<"[#6]-1(=[#6]-[#6](-c:2:c:c(:c(:n:c-1:2)-[#7](-[#1])-[#1])-[#6]#[#7])=[#6])-[#6]#[#7]">(mol);
    validate<"[#6]-1(=[#7]-[#7](-[#6](-[#16]-1)=[#6](-[#1])-[#6]:[#6])-[#6]:[#6])-[#6]=[#8]">(mol);
    validate<"[#6]-1(=[#8])-[#6](-[#6](-[#6]#[#7])=[#6](-[#1])-[#7])-[#6](-[#7])-[#6]=[#6]-1">(mol);
    validate<"[#6]-1(=[#8])-[#6](=[#6](-[#1])-[$([#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1),$([#6]:1:[#6]:[#6]:[#6]:[!#6&!#1]:1)])-[#7]=[#6](-[!#1]:[!#1]:[!#1])-[$([#16]),$([#7]-[!#1]:[!#1])]-1">(mol);
    validate<"[#6]-1-3=[#6](-[#6](-[#7]-c:2:c:c:c:c:c-1:2)(-[#6])-[#6])-[#16]-[#16]-[#6]-3=[!#1]">(mol);
}

template void Rarey_smarts_part_70<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_70<RDKit::ROMol>(RDKit::ROMol &mol);
