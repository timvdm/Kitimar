#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_68(Mol &mol)
{
    // SMARTS 336 - 340
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C-[$([CH]),$(C([CX4,c]))]=O">(mol);
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C-C#N">(mol);
    validate<"[OH][CX4][Cl,Br,I]">(mol);
    validate<"[Cl,Br,I][$([CX4][CH]=O),$([CX4]C(=O)[CX4,c])]">(mol);
    validate<"[OH][CX4][$([NH2]),$([NH][CX4]),$(N([CX4])[CX4])]">(mol);
}

template void RDKit_smarts_part_68<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_68<RDKit::ROMol>(RDKit::ROMol &mol);
