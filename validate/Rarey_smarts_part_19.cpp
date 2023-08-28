#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_19(Mol &mol)
{
    // SMARTS 91 - 95
    validate<"C=P">(mol);
    validate<"C=[OD2H0,SD2H0]-[SD1H0,OD1H0]">(mol);
    validate<"C=[SD1H0]">(mol);
    validate<"CC(C)=[CH][CH2][OH]">(mol);
    validate<"CC=C(C)[CH2][OH]">(mol);
}

template void Rarey_smarts_part_19<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_19<RDKit::ROMol>(RDKit::ROMol &mol);
