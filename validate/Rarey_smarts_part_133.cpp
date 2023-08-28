#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_133(Mol &mol)
{
    // SMARTS 661 - 665
    validate<"[C;$(C(=O))]-;!@[n]">(mol);
    validate<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;D3;!R;$(C(=O)([#6;R])!=[#6;!D1]~[!D1]);!$(C[O,S,N;!R])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol);
}

template void Rarey_smarts_part_133<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_133<RDKit::ROMol>(RDKit::ROMol &mol);
