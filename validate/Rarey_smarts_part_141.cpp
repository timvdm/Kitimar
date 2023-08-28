#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_141(Mol &mol)
{
    // SMARTS 701 - 705
    validate<"[CD2H1]=[ND2H0][ND2H1]">(mol);
    validate<"[CD2H1]=[ND2H0][OX2]">(mol);
    validate<"[CD2H1]=[OD1H0]">(mol);
    validate<"[CD2H2]-[CD3H0]=[ND1H1]">(mol);
    validate<"[CD3H0,R](=[SD1H0])([ND2H1,R])([ND2H1,R])">(mol);
}

template void Rarey_smarts_part_141<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_141<RDKit::ROMol>(RDKit::ROMol &mol);
