#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_246(Mol &mol)
{
    // SMARTS 1226 - 1230
    validate<"p-[F]">(mol);
    validate<"p-[I]">(mol);
    validate<"s1c(c(c-2c1-[#7](-[#1])-[#6](-[#6](=[#6]-2-[#1])-[#6](=[#8])-[#8]-[#1])=[#8])-[#7](-[#1])-[#1])-[#6](=[#8])-[#7]-[#1]">(mol);
    validate<"s1ccc(c1)-[#8]-[#1]">(mol);
    validate<"s1ccnc1-c2c(n(nc2-[#1])-[#1])-[#7](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_246<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_246<RDKit::ROMol>(RDKit::ROMol &mol);
