#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_153(Mol &mol)
{
    // SMARTS 761 - 765
    validate<"[CX4][CH]=[O,SX2]">(mol);
    validate<"[CX4v4]">(mol);
    validate<"[C](=O)([C,c,O,S])[C,c,O,S]">(mol);
    validate<"[Cl,Br,I][CH2]">(mol);
    validate<"[Cl]-c:2:c:c:1:n:o:n:c:1:c:c:2">(mol);
}

template void Rarey_smarts_part_153<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_153<RDKit::ROMol>(RDKit::ROMol &mol);
