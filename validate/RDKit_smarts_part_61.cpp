#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_61(Mol &mol)
{
    // SMARTS 301 - 305
    validate<"[S,s;D2]C[S,s;D2]">(mol);
    validate<"[$([S,s]~[S,s]~[C,c]=S),$([S,s]~[C,c](=S)~[S,s,N]),$([S;D2;R0]-S~O)]">(mol);
    validate<"[n+]-C">(mol);
    validate<"[NH]=C([NH2])c">(mol);
    validate<"O=CN[OH]">(mol);
}

template void RDKit_smarts_part_61<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_61<RDKit::ROMol>(RDKit::ROMol &mol);
