#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_192(Mol &mol)
{
    // SMARTS 956 - 960
    validate<"[r;!r3;!r4;!r5;!r6;!r7]">(mol);
    validate<"[sX2r5]">(mol);
    validate<"a-[Br]">(mol);
    validate<"a-[Cl]">(mol);
    validate<"a-[F]">(mol);
}

template void Rarey_smarts_part_192<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_192<RDKit::ROMol>(RDKit::ROMol &mol);
