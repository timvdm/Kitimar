#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_135(Mol &mol)
{
    // SMARTS 671 - 675
    validate<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol);
    validate<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]">(mol);
    validate<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]">(mol);
    validate<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;D3;R;!$([C]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol);
}

template void Rarey_smarts_part_135<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_135<RDKit::ROMol>(RDKit::ROMol &mol);
