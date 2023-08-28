#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_54(Mol &mol)
{
    // SMARTS 266 - 270
    validate<"[#6;X4]-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#6](=[#8])-[#7](-[#1])-[#6](-[#1])(-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#16]-[#6](-[#1])(-[#1])-[#1])-[#6](=[#8])-[#8]-[#1])-[#1])-[#1]">(mol);
    validate<"[#6X3](=[SX1])([!N])[!N]">(mol);
    validate<"[#6X3v3-1,#6X3v4+0,#6X3v3+1]">(mol);
    validate<"[#6]#[#6]-[#6](=[#8])-[#6]#[#6]">(mol);
    validate<"[#6](-[#16])(-[#7])=[#6](-[#1])-[#6]=[#6](-[#1])-[#6]=[#8]">(mol);
}

template void Rarey_smarts_part_54<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_54<RDKit::ROMol>(RDKit::ROMol &mol);
