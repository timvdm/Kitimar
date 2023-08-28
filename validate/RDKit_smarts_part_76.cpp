#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_76(Mol &mol)
{
    // SMARTS 376 - 380
    validate<"[$([#6]);!$(C=[O,S,N])]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol);
    validate<"[$(O([$([CX4,c]);!$(C[O,N,S])])[$([CX4,c]);!$(C[O,N,S])]);!$(O1[CX4][CX4]1)]">(mol);
    validate<"[OH][CH,CH2]O[CX4,c]">(mol);
    validate<"O([#6])-C([#6])([#6])-[OH]">(mol);
    validate<"[OH][$([NX3]([C;!$(C=[O,S,N])])[C;!$(C=[O,S,N])]),$([NH][CX4])]">(mol);
}

template void RDKit_smarts_part_76<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_76<RDKit::ROMol>(RDKit::ROMol &mol);
