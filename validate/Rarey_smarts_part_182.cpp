#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_182(Mol &mol)
{
    // SMARTS 906 - 910
    validate<"[SD4H0](=[OD1H0])(=[OD1H0])[Cl,Br,F]">(mol);
    validate<"[SH]">(mol);
    validate<"[SX2]">(mol);
    validate<"[SX4](C)(C)(=O)=N">(mol);
    validate<"[X4;R2;r4,r5,r6](@[r4,r5,r6])(@[r4,r5,r6])(@[r4,r5,r6])@[r4,r5,r6]">(mol);
}

template void Rarey_smarts_part_182<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_182<RDKit::ROMol>(RDKit::ROMol &mol);
