#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_75(Mol &mol)
{
    // SMARTS 371 - 375
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3][$([NH2]),$([NH][CX4]),$([N;R0]([CX4])[CX4])]);!$(CC=[O,S,N]);!$(C[O,S])]">(mol);
    validate<"[$([CX3]C(=O)[CX4,c]);!$(CC=[S,N]);!$(C[O,S,N])]=[$([CX3]C(=O)[CX4,c]);!$(CC=[S,N]);!$(C[N,S])]">(mol);
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3][OH]);!$(CC=[O,S,N]);!$(C[N,S])]">(mol);
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3]O[CX4]);!$(CC=[O,S,N]);!$(C[N,S])]">(mol);
    validate<"O1[CX4][CX4]1">(mol);
}

template void RDKit_smarts_part_75<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_75<RDKit::ROMol>(RDKit::ROMol &mol);
