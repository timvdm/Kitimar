#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_189(Mol &mol)
{
    // SMARTS 941 - 945
    validate<"[cR1]1[cR1][cR1][cR1][cR1][cR1]1">(mol);
    validate<"[nD3H0,R](~[OD1H0])(a)a">(mol);
    validate<"[nH0;!$(n-C);!$(n(:c)(:c):a)]1ccccc1">(mol);
    validate<"[nH1](n)nn">(mol);
    validate<"[nH1]1C(=O)CC(=O)O1">(mol);
}

template void Rarey_smarts_part_189<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_189<RDKit::ROMol>(RDKit::ROMol &mol);
