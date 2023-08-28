#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_108(Mol &mol)
{
    // SMARTS 536 - 540
    validate<"[#8]=[#6]-4-[#6]-[#6]-[#6]-3-[#6]-2-[#6](=[#8])-[#6]-[#6]-1-[#6]-[#6]-[#6]-[#6]-1-[#6]-2-[#6]-[#6]-[#6]-3=[#6]-4">(mol);
    validate<"[#8]=[#6]-[#6](-[#1])=[#6](-[#6]#[#7])-[#6]">(mol);
    validate<"[#8]=[#6]-[#6]-1=[#6](-[#16]-[#6](=[#6](-[#1])-[#6])-[#16]-1)-[#6]=[#8]">(mol);
    validate<"[#8]=[#6]-[#6]=[#6](-[#1])-[#8]-[#1]">(mol);
    validate<"[#8]=[#6]-[#7](-[#1])-c:1:c(-[#6]:[#6]):n:c(-[#6](-[#1])(-[#1])-[#6]#[#7]):s:1">(mol);
}

template void Rarey_smarts_part_108<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_108<RDKit::ROMol>(RDKit::ROMol &mol);
