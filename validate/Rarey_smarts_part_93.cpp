#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_93(Mol &mol)
{
    // SMARTS 461 - 465
    validate<"[#7](-[#1])-c1nc(nc2nnc(n12)-[#16]-[#6])-[#7](-[#1])-[#6]">(mol);
    validate<"[#7](-[#1])-c:1:n:c(:c:s:1)-c:2:c:n:c(-[#7](-[#1])-[#1]):s:2">(mol);
    validate<"[#7](-[#6](-[#1])-[#1])(-[#6](-[#1])-[#1])-[#6](-[#1])=[#7]-[#6](-[#6](-[#1])-[#1])=[#7]-[#7](-[#6](-[#1])-[#1])-[#6]:[#6]">(mol);
    validate<"[#7](-[#6](-[#1])-[#1])(-[#6](-[#1])-[#1])-c:1:c(:c(:c(:o:1)-[#6]=[#7]-[#7](-[#1])-[#6]=[!#6&!#1])-[#1])-[#1]">(mol);
    validate<"[#7](-c:1:c:c:c:c:c:1)-[#16](=[#8])(=[#8])-[#6]:2:[#6]:[#6]:[#6]:[#6]:3:[#7]:[$([#8]),$([#16])]:[#7]:[#6]:2:3">(mol);
}

template void Rarey_smarts_part_93<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_93<RDKit::ROMol>(RDKit::ROMol &mol);
