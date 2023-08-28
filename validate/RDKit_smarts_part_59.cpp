#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_59(Mol &mol)
{
    // SMARTS 291 - 295
    validate<"[Be,B,Al,Ti,Cr,Mn,Fe,Co,Ni,Cu,Pd,Ag,Sn,Pt,Au,Hg,Pb,Bi,As,Sb,Gd,Te]">(mol);
    validate<"O=CON1C(=O)CCC1=O">(mol);
    validate<"O=COn1cncc1">(mol);
    validate<"Fc1c(OC=O)c(F)c(F)c(F)c1F">(mol);
    validate<"[S;D2][C;R0](C)(C)C">(mol);
}

template void RDKit_smarts_part_59<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_59<RDKit::ROMol>(RDKit::ROMol &mol);
