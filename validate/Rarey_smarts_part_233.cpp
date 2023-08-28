#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_233(Mol &mol)
{
    // SMARTS 1161 - 1165
    validate<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7](-[#6;X4]-[#7]-3-[#1])-[#1]">(mol);
    validate<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7](-[#7]=[#6]-3)-[#1]">(mol);
    validate<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7]-[#6]=[#7]-3">(mol);
    validate<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7]-[#7]=[#7]-3">(mol);
    validate<"c:2:c:c:1:n:c(:c(:n:c:1:c:c:2)-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]:[#6])-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]:[#6]">(mol);
}

template void Rarey_smarts_part_233<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_233<RDKit::ROMol>(RDKit::ROMol &mol);
