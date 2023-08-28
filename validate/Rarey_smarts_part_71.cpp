#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_71(Mol &mol)
{
    // SMARTS 351 - 355
    validate<"[#6]-1-[#7](-[#1])-[#7](-[#1])-[#6](=[#16])-[#7]-[#7]-1-[#1]">(mol);
    validate<"[#6]-1:[#6]-[#6](=[#8])-[#6]=[#6]-1-[#7]=[#6](-[#1])-[#7](-[#6;X4])-[#6;X4]">(mol);
    validate<"[#6]-1:[#6]-[#7]=[#6]-[#6](=[#6]-[#7]-[#6])-[#16]-1">(mol);
    validate<"[#6]-1:[#6]-[#8]-[#6]-2-[#6](-[#1])(-[#1])-[#6](=[#8])-[#8]-[#6]-1-2">(mol);
    validate<"[#6]-1=,:[#6]-[#6](-[#6](-[$([#8]),$([#16])]-1)=[#6]-[#6]=[#8])=[#8]">(mol);
}

template void Rarey_smarts_part_71<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_71<RDKit::ROMol>(RDKit::ROMol &mol);
