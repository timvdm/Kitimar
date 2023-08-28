#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_73(Mol &mol)
{
    // SMARTS 361 - 365
    validate<"[#6]-2(-[#1])(-[#8]-[#1])-[#6]:1:[#7]:[!#6&!#1]:[#7]:[#6]:1-[#6](-[#1])(-[#8]-[#1])-[#6]=[#6]-2">(mol);
    validate<"[#6]-2(-[#6]=[#7]-c:1:c:c(:c:c:c:1-[#8]-2)-[Cl])=[#8]">(mol);
    validate<"[#6]-2(-[#6]=[#7]-c:1:c:c:c:c:c:1-[#7]-2)=[#6](-[#1])-[#6]=[#8]">(mol);
    validate<"[#6]-2(=[#16])-[#7](-[#6](-[#1])(-[#1])-c:1:c:c:c:o:1)-[#6](=[#7]-[#7]-2-[#1])-[#6]:[#6]">(mol);
    validate<"[#6]-2(=[#16])-[#7]-1-[#6]:[#6]-[#7]=[#7]-[#6]-1=[#7]-[#7]-2-[#1]">(mol);
}

template void Rarey_smarts_part_73<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_73<RDKit::ROMol>(RDKit::ROMol &mol);
