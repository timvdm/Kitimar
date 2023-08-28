#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_162(Mol &mol)
{
    // SMARTS 806 - 810
    validate<"[NH]([CX4])[CX4]">(mol);
    validate<"[NH](c)S(=O)=O">(mol);
    validate<"[NH]S(=O)(=O)C(F)(F)F">(mol);
    validate<"[NX1]#[CX2]">(mol);
    validate<"[NX2-]">(mol);
}

template void Rarey_smarts_part_162<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_162<RDKit::ROMol>(RDKit::ROMol &mol);
