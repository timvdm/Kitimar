#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_180(Mol &mol)
{
    // SMARTS 896 - 900
    validate<"[S-][CX3](=S)[#6]">(mol);
    validate<"[S;!H0]">(mol);
    validate<"[S;X4;$(S(=O)(=O))]-;!@[#7;X3;!D1;$([#7](S(=O)(=O))[#6;!D1][!D1]);!$([#7]C(=O))]">(mol);
    validate<"[S;X4;$(S(=O)(=O))]-;!@[N;X3;!$(N[N,O,S])]">(mol);
    validate<"[SD1H1]">(mol);
}

template void Rarey_smarts_part_180<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_180<RDKit::ROMol>(RDKit::ROMol &mol);
