#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_117(Mol &mol)
{
    // SMARTS 581 - 585
    validate<"[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N])]">(mol);
    validate<"[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]">(mol);
    validate<"[$([$([NX3]=O),$([NX3+][O-])])]">(mol);
    validate<"[$([$([NX4]=O),$([NX4+][O-])])]">(mol);
    validate<"[$([CH]=[CH2,CH]),$(C(C)=[CH2,CH]),$(C#C);!$(C(C)=CC)][CH2][OH]">(mol);
}

template void Rarey_smarts_part_117<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_117<RDKit::ROMol>(RDKit::ROMol &mol);
