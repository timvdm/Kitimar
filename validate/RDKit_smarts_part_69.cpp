#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_69(Mol &mol)
{
    // SMARTS 341 - 345
    validate<"O=C([c,CX4])[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol);
    validate<"[$([NH2][CX4]),$([$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]),$([ND3]([CX4])([CX4])[CX4])]">(mol);
    validate<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4]);!$(NC=O)][CX4]C(=O)[OH]">(mol);
    validate<"O(C(=O)[#6])(C(=O)[#6])">(mol);
    validate<"c[NH,NH2]">(mol);
}

template void RDKit_smarts_part_69<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_69<RDKit::ROMol>(RDKit::ROMol &mol);
