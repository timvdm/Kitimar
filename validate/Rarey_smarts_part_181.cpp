#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_181(Mol &mol)
{
    // SMARTS 901 - 905
    validate<"[SD2H0,OD2H0]-[SD2H0,OD2H0]">(mol);
    validate<"[SD2H0]-[Cl,Br,F]">(mol);
    validate<"[SD4H0](=[OD1H0])(=[OD1H0])(-[OD2H0]-C)(-[OD2H0]-C)">(mol);
    validate<"[SD4H0](=[OD1H0])(=[OD1H0])(-[OD2H0]-C)-C">(mol);
    validate<"[SD4H0](=[OD1H0])(=[OD1H0])(-a)-a">(mol);
}

template void Rarey_smarts_part_181<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_181<RDKit::ROMol>(RDKit::ROMol &mol);
