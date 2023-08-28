#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_47(Mol &mol)
{
    // SMARTS 231 - 235
    validate<"[#16;X2]-1-[#6]=[#6](-[#6]#[#7])-[#6](-[#6])(-[#6]=[#8])-[#6](=[#6]-1-[#7](-[#1])-[#1])-[$([#6]=[#8]),$([#6]#[#7])]">(mol);
    validate<"[#16X2H0]">(mol);
    validate<"[#16X2H0][!#16]">(mol);
    validate<"[#16X2H0][#16X2H0]">(mol);
    validate<"[#16X2H]">(mol);
}

template void Rarey_smarts_part_47<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_47<RDKit::ROMol>(RDKit::ROMol &mol);
