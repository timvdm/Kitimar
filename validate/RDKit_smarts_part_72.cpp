#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_72(Mol &mol)
{
    // SMARTS 356 - 360
    validate<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C[CX4]C(=O)[OH]">(mol);
    validate<"[CX4]([OH])[CX4]C(=O)[OH]">(mol);
    validate<"[CX4]([Cl,Br,I])[CX4]C(=O)[OH]">(mol);
    validate<"[#6]C(=O)[CX4]C(=O)[OH]">(mol);
    validate<"[#6]C(=O)[CX4]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol);
}

template void RDKit_smarts_part_72<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_72<RDKit::ROMol>(RDKit::ROMol &mol);
