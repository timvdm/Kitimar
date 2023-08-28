#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_175(Mol &mol)
{
    // SMARTS 871 - 875
    validate<"[OX2,OX1-][OX2,OX1-]">(mol);
    validate<"[OX2H+]=*">(mol);
    validate<"[OX2H,OX1H0-]">(mol);
    validate<"[OX2H0]">(mol);
    validate<"[OX2H]">(mol);
}

template void Rarey_smarts_part_175<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_175<RDKit::ROMol>(RDKit::ROMol &mol);
