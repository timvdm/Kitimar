#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_78(Mol &mol)
{
    // SMARTS 386 - 390
    validate<"O=C=N-[#6]">(mol);
    validate<"S=C=N-[#6]">(mol);
    validate<"O([CX4,c])-C([CX4,c])([CX4,c])-O([CX4,c])">(mol);
    validate<"[CX4,c]C(=O)[CX4,c]">(mol);
    validate<"[C;R0;X4]!@[CX4]!@[CX4]!@[CX4]!@[CX4]!@[CX4]!@[C;R0;X4]">(mol);
}

template void RDKit_smarts_part_78<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_78<RDKit::ROMol>(RDKit::ROMol &mol);
