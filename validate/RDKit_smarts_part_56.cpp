#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_56(Mol &mol)
{
    // SMARTS 276 - 280
    validate<"[N+]#N-*">(mol);
    validate<"[C,c]-[S;D2]-[O,N]">(mol);
    validate<"[Cl,Br,I,F][S,P,Si,N]">(mol);
    validate<"[Cl,Br,I]CC[S,N]">(mol);
    validate<"[Si]-O-*">(mol);
}

template void RDKit_smarts_part_56<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_56<RDKit::ROMol>(RDKit::ROMol &mol);
