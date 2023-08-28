#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_105(Mol &mol)
{
    // SMARTS 521 - 525
    validate<"[#8]-[#6](=[#8])-[#6](-[#1])(-[#1])-[#16;X2]-[#6](=[#7]-[#6]#[#7])-[#7](-[#1])-c:1:c:c:c:c:c:1">(mol);
    validate<"[#8]=[#16](=[#8])(-[#6]:[#6])-[#7](-[#1])-[#7](-[#1])-c1nc(cs1)-[#6]:[#6]">(mol);
    validate<"[#8]=[#16](=[#8])(-[#6]:[#6])-[#7](-[#1])-c1nc(cs1)-[#6]:[#6]">(mol);
    validate<"[#8]=[#16](=[#8])-[#6](-[#6]#[#7])=[#7]-[#7]-[#1]">(mol);
    validate<"[#8]=[#6](-c:1:c(:c(:n:c(:c:1-[#1])-[#8]-[#6](-[#1])(-[#1])-[#1])-[#8]-[#6](-[#1])(-[#1])-[#1])-[#1])-[#7](-[#1])-[#6](-[#1])(-[#6](-[#1])-[#1])-[#6](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_105<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_105<RDKit::ROMol>(RDKit::ROMol &mol);
