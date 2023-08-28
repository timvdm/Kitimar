#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_124(Mol &mol)
{
    // SMARTS 616 - 620
    validate<"[$([SX4](=O)(=O)(O)O),$([SX4+2]([O-])([O-])(O)O)]">(mol);
    validate<"[$([SX4](=[OX1])(=[OX1])([!O])[NX3]),$([SX4+2]([OX1-])([OX1-])([!O])[NX3])]">(mol);
    validate<"[$([cX2+](:*):*)]">(mol);
    validate<"[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*)]">(mol);
    validate<"[$([cX3](:*):*),$([cX2+](:*):*)]">(mol);
}

template void Rarey_smarts_part_124<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_124<RDKit::ROMol>(RDKit::ROMol &mol);
