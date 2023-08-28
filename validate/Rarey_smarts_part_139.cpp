#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_139(Mol &mol)
{
    // SMARTS 691 - 695
    validate<"[C;X4;!D4]-;!@[C;D3]">(mol);
    validate<"[C;X4;!D4]-;!@[C;X4;!D4;!R]">(mol);
    validate<"[C;X4;!D4]-;!@[N;X3;!D3]">(mol);
    validate<"[C;X4;!D4]-;!@[O;D2;$(O[C,c]);!$(OCO);!$(OC=*)]">(mol);
    validate<"[C;X4;!D4]-;!@[c;$(c[n&H0,cH,S,O])]">(mol);
}

template void Rarey_smarts_part_139<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_139<RDKit::ROMol>(RDKit::ROMol &mol);
