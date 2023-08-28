#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_186(Mol &mol)
{
    // SMARTS 926 - 930
    validate<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[c;X3;R;!$([c]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol);
    validate<"[c;X3;R;$([c]@[N,n,S,s,O,o])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([#6;!D1;R])[#6;!D1;R]!=,@[!D1])]">(mol);
}

template void Rarey_smarts_part_186<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_186<RDKit::ROMol>(RDKit::ROMol &mol);
