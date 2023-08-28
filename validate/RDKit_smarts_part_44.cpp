#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_44(Mol &mol)
{
    // SMARTS 216 - 220
    validate<"N1-C-O1">(mol);
    validate<"[C;H2]-O-S(=O)(=O)-C">(mol);
    validate<"C1-C-O1">(mol);
    validate<"c1-[N;a]-c([F,Br,I,Cl])-c-c-c1">(mol);
    validate<"C=[C;H0]([F,Br,I,Cl])([F,Br,I,Cl])">(mol);
}

template void RDKit_smarts_part_44<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_44<RDKit::ROMol>(RDKit::ROMol &mol);
