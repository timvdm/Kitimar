#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_131(Mol &mol)
{
    // SMARTS 651 - 655
    validate<"[C;!D4;!D1;!R;$(C(=O));$(C([#7;R;D3])[!S;!D1][!D1])]-;!@[#7;R;D3;$([#7]@[#6](=O))]">(mol);
    validate<"[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[!N;!D1][!D1])]-;!@[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!N;!D1][!D1])]">(mol);
    validate<"[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[!S;!D1][!D1])]-;!@[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!S;!D1][!D1])]">(mol);
    validate<"[C;!D4;!D1;!R;$(C(=O));$(C([O;D2])[#6;!D1]!=[!D1])]-;!@[O;D2;$(O(C(=O))[#6;!D1][#6;!D1]);!$(OCO);!$(O[P,S])]">(mol);
    validate<"[C;!D4;!D1;!R;$(C(=O));$(C([n])[!S;!D1][!D1])]-;!@[n;!$(n@[#6](=O))]">(mol);
}

template void Rarey_smarts_part_131<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_131<RDKit::ROMol>(RDKit::ROMol &mol);
