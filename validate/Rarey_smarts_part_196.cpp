#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_196(Mol &mol)
{
    // SMARTS 976 - 980
    validate<"c-@P">(mol);
    validate<"c-@S">(mol);
    validate<"c1(O[CH3])ccc([CH2][CH]=[CH2])cc1">(mol);
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])c([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc1">(mol);
    validate<"c1([F,Cl,Br,I,$(N(=O)~O)])c([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cccc1([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])">(mol);
}

template void Rarey_smarts_part_196<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_196<RDKit::ROMol>(RDKit::ROMol &mol);
