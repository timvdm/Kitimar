#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_203(Mol &mol)
{
    // SMARTS 1011 - 1015
    validate<"c2(nc:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1])n2-[#6])-[#7](-[#1])-[#6](-[#7](-[#1])-c:3:c(:c:c:c:c:3-[#1])-[#1])=[#8]">(mol);
    validate<"c2-3:c:c:c:1:c:c:c:c:c:1:c2-[#6](-[#1])-[#6;X4]-[#7]-[#6]-3=[#6](-[#1])-[#6](=[#8])-[#7](-[#6](-[#1])-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"c3cn1c(nc(c1-[#7]-[#6])-c:2:c:c:c:c:n:2)cc3">(mol);
    validate<"c:1(:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#6](=[#6](-[#1])-[#1])-[#6](-[#1])-[#1])-[#1])-[#6](-[#6;X4])(-[#6;X4])-[#7](-[#1])-[#6](=[#8])-[#7](-[#6](-[#1])(-[#1])-[#6](-[#1])-[#1])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#6](-[#1])-[#6](-[#1])(-[#1])-[#6]:[#6]">(mol);
    validate<"c:1(:c(:c(:c(:c(:c:1-[#1])-[#1])-[#6](-[#1])-[#1])-[#7](-[#1])-[#6](=[#8])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#6]:[#6])-[#1])-[#7](-[#1])-[#6](=[#8])-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#6]:[#6]">(mol);
}

template void Rarey_smarts_part_203<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_203<RDKit::ROMol>(RDKit::ROMol &mol);
