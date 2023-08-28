#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_200(Mol &mol)
{
    // SMARTS 996 - 1000
    validate<"c1c(=[O,NH2,NH])c(=[O,NH2,NH])ccc1">(mol);
    validate<"c1c(=[O,NH2,NH])ccc(=[O,NH2,NH])c1">(mol);
    validate<"c1c([OH,NH2,NH])c([OH,NH2,NH,$(N=N),$(N(C)C)])ccc1">(mol);
    validate<"c1c([OH,NH2,NH])cc([OH,NH2,NH,$(N(C)C)])cc1">(mol);
    validate<"c1c([OH,NH2,NH])ccc([OH,NH2,NH,$(N=N),$(N(C)C)])c1">(mol);
}

template void Rarey_smarts_part_200<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_200<RDKit::ROMol>(RDKit::ROMol &mol);
