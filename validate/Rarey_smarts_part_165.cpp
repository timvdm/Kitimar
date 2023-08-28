#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_165(Mol &mol)
{
    // SMARTS 821 - 825
    validate<"[NX3]([CX4])([CX4])[CX4]">(mol);
    validate<"[NX3]-[Cl,Br,F,I]">(mol);
    validate<"[NX3]-[OX2]">(mol);
    validate<"[NX3][$(C=C),$(cc)]">(mol);
    validate<"[NX3][CX2]#[NX1]">(mol);
}

template void Rarey_smarts_part_165<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_165<RDKit::ROMol>(RDKit::ROMol &mol);
