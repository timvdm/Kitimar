#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_184(Mol &mol)
{
    // SMARTS 916 - 920
    validate<"[c;$(c([n&H0,cH,S,O]))]-;!@[N;X3;!$(N[N,O,S])]">(mol);
    validate<"[c;$(c[n&H0,cH,S,O])]-;!@[c;$(c[n&H0,cH,S,O])]">(mol);
    validate<"[c;$(c[n&H0,cH,S,O])]-;!@[n]">(mol);
    validate<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol);
    validate<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]">(mol);
}

template void Rarey_smarts_part_184<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_184<RDKit::ROMol>(RDKit::ROMol &mol);
