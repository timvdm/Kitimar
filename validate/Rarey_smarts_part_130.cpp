#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_130(Mol &mol)
{
    // SMARTS 646 - 650
    validate<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([n])!=[#6;!D1]~[!D1])]-;!@[n;!$(n@[#6](=O))]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);$(C([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol);
}

template void Rarey_smarts_part_130<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_130<RDKit::ROMol>(RDKit::ROMol &mol);
