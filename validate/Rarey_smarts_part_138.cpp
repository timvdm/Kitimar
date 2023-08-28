#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_138(Mol &mol)
{
    // SMARTS 686 - 690
    validate<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[c;X3;R;$([c]@[N,n,S,s,O,o])]">(mol);
    validate<"[C;D3;R;$([C]@[N,n,S,s,O,o])]-;!@[n;!$(n@[#6](=O))]">(mol);
    validate<"[C;H1](=[O,S])[C,c]">(mol);
    validate<"[C;H2;$([CH2]([CH2])[CH2])][H]">(mol);
    validate<"[C;H3;$(Cc([cH])([cH,n&H0,O]))][H]">(mol);
}

template void Rarey_smarts_part_138<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_138<RDKit::ROMol>(RDKit::ROMol &mol);
