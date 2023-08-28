#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_63(Mol &mol)
{
    // SMARTS 311 - 315
    validate<"[$(N#CSc1sc(nc1)N),$([S,Se]1C(N)C(=O)[#6][#6]1)]">(mol);
    validate<"[!$(O=[C,S])][N;R0]=[C;R0]([C,c])[C,c]">(mol);
    validate<"O=CC([Cl,Br,I,F])([Cl,Br,I,F])[Cl,Br,I,F]">(mol);
    validate<"[#6]SC(=O)">(mol);
    validate<"[$(O([CX4,c])!@[CH,CH2]!@O[CX4,c])]">(mol);
}

template void RDKit_smarts_part_63<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_63<RDKit::ROMol>(RDKit::ROMol &mol);
