#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_176(Mol &mol)
{
    // SMARTS 876 - 880
    validate<"[OX2H]P">(mol);
    validate<"[OX2H][#6X3]=[#6]">(mol);
    validate<"[OX2H][$(C=C),$(cc)]">(mol);
    validate<"[OX2H][CX3]=[OX1]">(mol);
    validate<"[OX2H][cX3]:[c]">(mol);
}

template void Rarey_smarts_part_176<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_176<RDKit::ROMol>(RDKit::ROMol &mol);
