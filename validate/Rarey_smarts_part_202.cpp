#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_202(Mol &mol)
{
    // SMARTS 1006 - 1010
    validate<"c1csc(n1)-[#7]-[#7]-[#16](=[#8])=[#8]">(mol);
    validate<"c1nc([F,Cl,Br,I,$(N(=O)~O)])ncn1">(mol);
    validate<"c2(c(-[#1])n(-[#6](-[#1])-[#1])c:3:c(:c(:c:1n(c(c(c:1:c2:3)-[#1])-[#1])-[#6](-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])-[#8]-[#6](-[#1])-[#1])-[#1]">(mol);
    validate<"c2(c(-[#7](-[#1])-[#1])n(-c:1:c:c:c:c:c:1-[#6](=[#8])-[#8]-[#1])nc2-[#6]=[#8])-[$([#6]#[#7]),$([#6]=[#16])]">(mol);
    validate<"c2(c-1n(-[#6](-[#6]=[#6]-[#7]-1)=[#8])nc2-c3cccn3)-[#6]#[#7]">(mol);
}

template void Rarey_smarts_part_202<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_202<RDKit::ROMol>(RDKit::ROMol &mol);
