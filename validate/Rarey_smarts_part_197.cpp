#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_197(Mol &mol)
{
    // SMARTS 981 - 985
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])ncc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc1">(mol);
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])ncc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])nc1">(mol);
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])ncccc1([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])">(mol);
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])ncccn1">(mol);
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])ncncc1">(mol);
}

template void Rarey_smarts_part_197<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_197<RDKit::ROMol>(RDKit::ROMol &mol);
