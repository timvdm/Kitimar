#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_140(Mol &mol)
{
    // SMARTS 696 - 700
    validate<"[C;X4;!D4]-;!@[n]">(mol);
    validate<"[CD1H2,CD2H1,CD3H0;!$(C=,#*)]">(mol);
    validate<"[CD1H2]=[CD2H1]-[OD2H0,SD2H0]-C">(mol);
    validate<"[CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0]">(mol);
    validate<"[CD2H1]=[CD3H0]-[ND1H2]">(mol);
}

template void Rarey_smarts_part_140<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_140<RDKit::ROMol>(RDKit::ROMol &mol);
