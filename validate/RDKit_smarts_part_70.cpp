#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_70(Mol &mol)
{
    // SMARTS 346 - 350
    validate<"c-[ND3]([#6])[#6]">(mol);
    validate<"c">(mol);
    validate<"c-[F,Cl,Br,I]">(mol);
    validate<"c-[Br,I]">(mol);
    validate<"Fcaa[F,Cl,Br,I,$([C,N,S]=O)]">(mol);
}

template void RDKit_smarts_part_70<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_70<RDKit::ROMol>(RDKit::ROMol &mol);
