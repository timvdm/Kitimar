#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_150(Mol &mol)
{
    // SMARTS 746 - 750
    validate<"[CX3](=O)[O-]">(mol);
    validate<"[CX3](=O)[OX1H0-,OX2H1]">(mol);
    validate<"[CX3](=O)[OX2H1]">(mol);
    validate<"[CX3](=[OX1])(O)O">(mol);
    validate<"[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]">(mol);
}

template void Rarey_smarts_part_150<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_150<RDKit::ROMol>(RDKit::ROMol &mol);
