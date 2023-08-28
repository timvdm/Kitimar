#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_226(Mol &mol)
{
    // SMARTS 1126 - 1130
    validate<"c:1:c:c(:c:c:c:1-[#7](-[#1])-[#1])-[#7](-[#6;X3])-[#6;X3]">(mol);
    validate<"c:1:c:c(:c:c:c:1-[#7](-[#6;X4])-[#6;X4])-[#6;X4]-[$([#8]-[#1]),$([#6]=[#6]-[#1]),$([#7]-[#6;X4])]">(mol);
    validate<"c:1:c:c(:c:c:c:1-[#7](-[#6;X4])-[#6;X4])-[#6]=[#6]">(mol);
    validate<"c:1:c:c(:c:c:c:1-[#8]-[#1])-[#7](-[#1])-[#16](=[#8])=[#8]">(mol);
    validate<"c:1:c:c(:c:c:c:1-[#8]-[#6;X4])-[#7](-[#6;X4])-[$([#1]),$([#6;X4])]">(mol);
}

template void Rarey_smarts_part_226<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_226<RDKit::ROMol>(RDKit::ROMol &mol);
