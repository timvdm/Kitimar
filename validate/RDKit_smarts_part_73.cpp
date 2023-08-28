#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_73(Mol &mol)
{
    // SMARTS 361 - 365
    validate<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol);
    validate<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]C(=O)[OH]">(mol);
    validate<"C(=O)[OH]">(mol);
    validate<"O([CX4,c])C(=O)O[CX4,c]">(mol);
    validate<"[$([CX4,c][CH]=O),$([CX4,c]C(=O)[CX4,c])]">(mol);
}

template void RDKit_smarts_part_73<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_73<RDKit::ROMol>(RDKit::ROMol &mol);
