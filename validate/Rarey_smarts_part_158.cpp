#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_158(Mol &mol)
{
    // SMARTS 786 - 790
    validate<"[N;X4]-;!@[C;D3;R;$([C]@[N,n,S,s,O,o])]">(mol);
    validate<"[N;X4]-;!@[c;X3;R;!$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[N;X4]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[ND1H2]a">(mol);
    validate<"[ND1H2]aaaa[ND1H2]">(mol);
}

template void Rarey_smarts_part_158<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_158<RDKit::ROMol>(RDKit::ROMol &mol);
