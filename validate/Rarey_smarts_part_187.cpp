#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_187(Mol &mol)
{
    // SMARTS 931 - 935
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1])[#6;!D1;R]!=,@[!D1])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1;R][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[O;D2;$(O([#6;!D1;R])[#6;!D1][!D1]);!$(OCO);!$(OC=;!@*);!$(O[P,S])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[S;$(S(=O)([#6;!D1;R])[#6;!D1;R][!D1])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[S;D2;$(S([#6;!D1;R])[#6;!D1;R][!D1]);!$(SCS);!$(SC=;!@*);!$(S[P,S])]">(mol);
}

template void Rarey_smarts_part_187<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_187<RDKit::ROMol>(RDKit::ROMol &mol);
