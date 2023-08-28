#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_51(Mol &mol)
{
    // SMARTS 251 - 255
    validate<"[#16]=[#6]-2-[#7](-[#1])-[#7]=[#6](-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])-[#1])-[#1])-[#8]-2">(mol);
    validate<"[#16]=[#6]-[#6](-[#6](-[#1])-[#1])=[#6](-[#6](-[#1])-[#1])-[#7](-[#6](-[#1])-[#1])-[#6](-[#1])-[#1]">(mol);
    validate<"[#16]=[#6]-c:1:c:c:c:2:c:c:c:c:n:1:2">(mol);
    validate<"[#6!H0,#1]">(mol);
    validate<"[#6!H0]~@[#6!H0]~@[#6!H0]~@[#6!H0]">(mol);
}

template void Rarey_smarts_part_51<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_51<RDKit::ROMol>(RDKit::ROMol &mol);
