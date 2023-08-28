#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_83(Mol &mol)
{
    // SMARTS 411 - 415
    validate<"[#6]=[#6]-[#6](-[#6]#[#7])(-[#6]#[#7])-[#6](-[#6]#[#7])=[#6]-[#7](-[#1])-[#1]">(mol);
    validate<"[#6]=[#6]-[#6](=[#8])-[#7]-c:1:c(:c(:c(:s:1)-[#6](=[#8])-[#8])-[#6]-[#1])-[#6]#[#7]">(mol);
    validate<"[#6]=[#7;!R]-c:1:c:c:c:c:c:1-[#8]-[#1]">(mol);
    validate<"[#6]C(=[O,SX2])C(=[O,SX2])[#6]">(mol);
    validate<"[#6]C(=[O,SX2])[CX4]C(=[O,SX2])[#6]">(mol);
}

template void Rarey_smarts_part_83<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_83<RDKit::ROMol>(RDKit::ROMol &mol);
