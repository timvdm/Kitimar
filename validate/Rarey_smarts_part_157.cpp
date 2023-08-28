#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_157(Mol &mol)
{
    // SMARTS 781 - 785
    validate<"[N;R0]=[N;R0]C#N">(mol);
    validate<"[N;R0]=[N;R0]CC=O">(mol);
    validate<"[N;R0][N;R0]C(=O)">(mol);
    validate<"[N;X4]-;!@[C;!D4;!D1;!R;!$(C(=O));$(C([N;X4;!D1])[#6;!D1][!D1])]">(mol);
    validate<"[N;X4]-;!@[C;D3;R;!$([C]@[N,n,S,s,O,o])]">(mol);
}

template void Rarey_smarts_part_157<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_157<RDKit::ROMol>(RDKit::ROMol &mol);
