#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_79(Mol &mol)
{
    // SMARTS 391 - 395
    validate<"[$([C;R1]);!$(C(N)N)](=O)@[$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol);
    validate<"[$([C;R1]);!$(C(O)N);!$(C(O)O)](=O)@[$(O);!$(O(C=O))]">(mol);
    validate<"OC(=O)CC(=O)[OH]">(mol);
    validate<"[r8,r9,r10,r11,r12,r13,r14]">(mol);
    validate<"*-C#N">(mol);
}

template void RDKit_smarts_part_79<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_79<RDKit::ROMol>(RDKit::ROMol &mol);
