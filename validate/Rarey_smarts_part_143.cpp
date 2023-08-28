#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_143(Mol &mol)
{
    // SMARTS 711 - 715
    validate<"[CD3H0]=[CD2H0]=[OD1H0]">(mol);
    validate<"[CD3H1;r6;R2,R3]12-aa-[CD3H1;r6;R2,R3](-aa1)aa2">(mol);
    validate<"[CD3H1]([CD3H0](~[OD1])(~[OD1]))([ND1H2])[CD2H2](a1aaaa2aaaaa12)">(mol);
    validate<"[CH,CH2,CH3;!$([CH2]CC=[O,S])][F,Cl,Br,I,$(OS(=O)(=O)[#6,#1]),$(OS(=O)(=O)O[#6,#1])]">(mol);
    validate<"[CH2,CH3][NX3][NX2]=[O,S]">(mol);
}

template void Rarey_smarts_part_143<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_143<RDKit::ROMol>(RDKit::ROMol &mol);
