#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_129(Mol &mol)
{
    // SMARTS 641 - 645
    validate<"[C,c](=N)N[C,S](=O)">(mol);
    validate<"[C;!D4;!D1;!R;!$(C(=O));$(C([N;X3;!D1])[#6;!D1][!D1])]-;!@[N;X3;!D1;!$(NC=;!@[N,O,S]);$(N([C;!D1;!R])[#6;!D1]!=[!D1])]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C(=O));$(C([O;D2])[#6;!D1]~[!D1])]-;!@[O;D2;$(O([C;!D1;!R])[#6;!D1][!D1]);!$(OCO);!$(OC=*);!$(O[P,S])]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C(=S));$(C([S;D2])[#6;!D1]~[!D1])]-;!@[S;D2;$(S([C;!D1;!R])[#6;!D1][!D1]);!$(SCS);!$(SC=*);!$(S[P,S])]">(mol);
    validate<"[C;!D4;!D1;!R;!$(C=[C,O,S,N]);!$(C[O,S,N;!R]);$(C([#7;R;D3])!=[#6;!D1]~[!D1])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol);
}

template void Rarey_smarts_part_129<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_129<RDKit::ROMol>(RDKit::ROMol &mol);
