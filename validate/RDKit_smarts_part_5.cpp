#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_5(Mol &mol)
{
    // SMARTS 21 - 25
    validate<"[N&D2](=O)">(mol);
    validate<"[P,S][Cl,Br,F,I]">(mol);
    validate<"N=C=N">(mol);
    validate<"[N+]#[C-]">(mol);
    validate<"C(=O)N(C(=O))OC(=O)">(mol);
}

template void RDKit_smarts_part_5<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_5<RDKit::ROMol>(RDKit::ROMol &mol);
