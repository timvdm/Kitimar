#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_54(Mol &mol)
{
    // SMARTS 266 - 270
    validate<"S(=O)(=O)-[C;H2]([F,Br,I,Cl])">(mol);
    validate<"[$(C-[Cl,Br,I]),$(O=C-[Cl,Br,I,F]),$(O=C([CH,CH2][Cl,Br,I,F])[O,C]),$(C~O~[Cl,Br,I,F][CH,CH2]),$(n1c([Cl,Br,I,F])nccc1);!$(C=C-[Cl,Br,I]);!$(ClC-[Cl,Br,I,F])]">(mol);
    validate<"[O,N,S]1CC1">(mol);
    validate<"O=COC=O">(mol);
    validate<"N=C=O">(mol);
}

template void RDKit_smarts_part_54<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_54<RDKit::ROMol>(RDKit::ROMol &mol);
