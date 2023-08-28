#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_188(Mol &mol)
{
    // SMARTS 936 - 940
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol);
    validate<"[cD3H0;r5,r6;R2](:[cD3H0;r5,r6;R2])(:[cD3H0;r5,r6;R2])">(mol);
    validate<"[cD3H0;r6;R2][*;r5,r6;R1][cD3H0;r6;R2]">(mol);
}

template void Rarey_smarts_part_188<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_188<RDKit::ROMol>(RDKit::ROMol &mol);
