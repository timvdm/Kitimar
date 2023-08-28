#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_198(Mol &mol)
{
    // SMARTS 986 - 990
    validate<"c1([OH])c(O[CH3])cccc1">(mol);
    validate<"c1([OH])ccc(O[CH3])cc1">(mol);
    validate<"c1([OH])ccc([CH2][CH]=[CH2])cc1">(mol);
    validate<"c1(c(c(c(n1-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#1])-[#1])-[#1])-[#6](-[#1])-[#1])-[#1])-[#6](=[#8])-[#8]-[#1]">(mol);
    validate<"c1(c-2c(c(n1-[#6](-[#8])=[#8])-[#6](-[#1])-[#1])-[#16]-[#6](-[#1])(-[#1])-[#16]-2)-[#6](-[#1])-[#1]">(mol);
}

template void Rarey_smarts_part_198<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_198<RDKit::ROMol>(RDKit::ROMol &mol);
