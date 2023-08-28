#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_199(Mol &mol)
{
    // SMARTS 991 - 995
    validate<"c1(coc(c1-[#1])-[#6](=[#16])-[#7]-2-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[!#1]-[#6](-[#1])(-[#1])-[#6]-2(-[#1])-[#1])-[#1]">(mol);
    validate<"c1(nn(c(c1-[$([#1]),$([#6]-[#1])])-[#8]-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#1])-[#1])-[#1])-[#6;X4]">(mol);
    validate<"c12ccccc1cccc2">(mol);
    validate<"c12ccccc1ccn2">(mol);
    validate<"c1c(-[#7](-[#1])-[#1])nnc1-c2c(-[#6](-[#1])-[#1])oc(c2-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_199<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_199<RDKit::ROMol>(RDKit::ROMol &mol);
