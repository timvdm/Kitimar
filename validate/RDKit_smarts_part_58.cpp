#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_58(Mol &mol)
{
    // SMARTS 286 - 290
    validate<"N([H,C;X4])([H,C;X4])-[C;R0;X4]-N[C;X4]([H,C;X4])([H,C;X4])">(mol);
    validate<"N#C[C;R0;X4]O[!$(O=[C,S])]">(mol);
    validate<"[C;R0](=[C;R0])-[S,O,N;R0][!$(O=[C,S])]">(mol);
    validate<"[O-]-[C;R0]=[C;R0]">(mol);
    validate<"[C;R0](=[C;R0])-[OH]">(mol);
}

template void RDKit_smarts_part_58<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_58<RDKit::ROMol>(RDKit::ROMol &mol);
