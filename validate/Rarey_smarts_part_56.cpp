#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_56(Mol &mol)
{
    // SMARTS 276 - 280
    validate<"[#6](-[#1])(-[#1])-[#16;X2]-c:1:n:c(:c(:n:1-!@[#6](-[#1])-[#1])-c:2:c:c:c:c:c:2)-[#1]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#16;X2]-c:1:n:n:c(:c(:n:1)-c:2:c(:c(:c(:o:2)-[#1])-[#1])-[#1])-c:3:c(:c(:c(:o:3)-[#1])-[#1])-[#1]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#16;X2]-c:2:n:n:c:1-[#6]:[#6]-[#7]=[#6]-[#8]-c:1:n:2">(mol);
    validate<"[#6](-[#1])(-[#1])-[#16]-[#6](=[#16])-[#7](-[#1])-[#6](-[#1])(-[#1])-[#6]:[#6]">(mol);
    validate<"[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#16]-[#6](-[#1])(-[#1])-c1cn(cn1)-[#1]">(mol);
}

template void Rarey_smarts_part_56<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_56<RDKit::ROMol>(RDKit::ROMol &mol);
