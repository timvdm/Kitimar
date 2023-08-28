#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_132(Mol &mol)
{
    // SMARTS 656 - 660
    validate<"[C;!D4;!D1;!R;$(C(S=O)[#6;!D1][!D1])]-;!@[S;$(S(=O)([C;!D1;!R])[#6;!D1][!D1])]">(mol);
    validate<"[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]~[C;!R]">(mol);
    validate<"[C;!r5]([C;!r5])=[C;!r5](C)[C;!r5]=[O,SX2;!r5]">(mol);
    validate<"[C;$(C(=O))]-;!@[N;!$(N(C=O)(C=O));!$(N[N,O,S])]">(mol);
    validate<"[C;$(C(=O))]-;!@[O;D2;$(O[C,c]);!$(OCO);!$(OC=*)]">(mol);
}

template void Rarey_smarts_part_132<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_132<RDKit::ROMol>(RDKit::ROMol &mol);
