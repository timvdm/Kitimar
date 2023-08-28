#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_76(Mol &mol)
{
    // SMARTS 376 - 380
    validate<"[#6]-[#6](=[#8])-[#6](-[#1])(-[#1])-[#16;X2]-c:3:n:n:c:2:c:1:c(:c(:c(:c(:c:1:n(:c:2:n:3)-[#1])-[#1])-[#1])-[#1])-[#1]">(mol);
    validate<"[#6]-[#6](=[#8])-[#6](-[#1])=[#6](-[#7](-[#1])-[#6])-[#6](=[#8])-[#8]-[#6]">(mol);
    validate<"[#6]-[#6]=[#6](-[F,Cl,Br,I])-[#6](=[#8])-[#6]">(mol);
    validate<"[#6]-[#7]-1-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#7](-[#6](-[#1])(-[#1])-[#6]-1(-[#1])-[#1])-[#7]=[#6](-[#1])-[#6]:[!#1]">(mol);
    validate<"[#6]1=,:[#6]C(=[O,SX1])N[SX2]1">(mol);
}

template void Rarey_smarts_part_76<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_76<RDKit::ROMol>(RDKit::ROMol &mol);
