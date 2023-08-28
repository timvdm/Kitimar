#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_80(Mol &mol)
{
    // SMARTS 396 - 400
    validate<"[#6]:[#6]-[#6](=[#16;X1])-[#16;X2]-[#6](-[#1])-[$([#6](-[#1])-[#1]),$([#6]:[#6])]">(mol);
    validate<"[#6]:[#6]-[#6](=[#8])-[#7](-[#1])-[#6](=[#8])-[#6](-[#6]#[#7])=[#6](-[#1])-[#7](-[#1])-[#6]:[#6]">(mol);
    validate<"[#6]:[#6]-[#7;!R]=[#6]-2-[#6](=[!#6&!#1])-c:1:c:c:c:c:c:1-[#7]-2">(mol);
    validate<"[#6]:[#6]-[#7](-[#1])-[#16](=[#8])(=[#8])-[#7](-[#1])-[#6]:[#6]">(mol);
    validate<"[#6]:[#6]-[#7](-[#1])-[#6](=[#8])-c1c(snn1)-[#7](-[#1])-[#6]:[#6]">(mol);
}

template void Rarey_smarts_part_80<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_80<RDKit::ROMol>(RDKit::ROMol &mol);
