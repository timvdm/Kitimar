#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_201(Mol &mol)
{
    // SMARTS 1001 - 1005
    validate<"c1c([OH])ccc(C=C[CH3])c1">(mol);
    validate<"c1ccc2c(c1)ccc3ccccc23">(mol);
    validate<"c1ccc2cc3ccccc3cc2c1">(mol);
    validate<"c1cncc1">(mol);
    validate<"c1csc(c1-[#7](-[#1])-[#1])-[#6](-[#1])=[#6](-[#1])-c2cccs2">(mol);
}

template void Rarey_smarts_part_201<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_201<RDKit::ROMol>(RDKit::ROMol &mol);
