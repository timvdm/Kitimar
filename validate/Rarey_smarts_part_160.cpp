#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_160(Mol &mol)
{
    // SMARTS 796 - 800
    validate<"[ND2H0][OD1H0]">(mol);
    validate<"[ND2H1,ND3H0](a1aaa2aaaaa2a1)a1aaaaa1">(mol);
    validate<"[ND3H0](~[OD1H0])(~[OD1H0])-A">(mol);
    validate<"[ND3H0](~[OD1H0])(~[OD1H0])-a">(mol);
    validate<"[NH,O,S;r3]">(mol);
}

template void Rarey_smarts_part_160<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_160<RDKit::ROMol>(RDKit::ROMol &mol);
