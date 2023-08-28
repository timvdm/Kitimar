#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_156(Mol &mol)
{
    // SMARTS 776 - 780
    validate<"[N+]#[C-]">(mol);
    validate<"[N,n;+1]">(mol);
    validate<"[N;!D1;!$(N(C=O)(C=O));$(N(C(=O))[!D1][!D1])]-;!@[C;!D4;!D1;!R;$(C(=O));$(C([N;!D1])[N;!D1][!D1])]">(mol);
    validate<"[N;!H0;!$([N;!H0]#C);!$([NH](C(F)(F)(F))S(=O)=O);!$([nH]1nnnc1);!$([nH]1nncn1);!$([nH]1ncnn1)]">(mol);
    validate<"[N;+0,+1;$(N(=O)~[O;H0;-0,-1])]">(mol);
}

template void Rarey_smarts_part_156<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_156<RDKit::ROMol>(RDKit::ROMol &mol);
