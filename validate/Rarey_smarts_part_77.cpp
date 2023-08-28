#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_77(Mol &mol)
{
    // SMARTS 381 - 385
    validate<"[#6]1C(=[O,S])[O,N,S][#6]1">(mol);
    validate<"[#6]1[#6](=[N,O,S])[#7,#8,#16][#6][#7,#8,#16]1">(mol);
    validate<"[#6]1[O,N,SX2][#6]1">(mol);
    validate<"[#6]:1(:[#7]:[#6](:[#7]:[!#1]:[#7]:1)-c:2:c(:c(:c(:o:2)-[#1])-[#1])-[#1])-[#16]-[#6;X4]">(mol);
    validate<"[#6]:1:2:[!#1]:[#7+](:[!#1]:[#6](:[!#1]:1:[#6]:[#6]:[#6]:[#6]:2)-[*])~[#6]:[#6]">(mol);
}

template void Rarey_smarts_part_77<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_77<RDKit::ROMol>(RDKit::ROMol &mol);
