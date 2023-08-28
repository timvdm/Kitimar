#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_123(Mol &mol)
{
    // SMARTS 611 - 615
    validate<"[$([OX2])]">(mol);
    validate<"[$([OX2]C=N)]">(mol);
    validate<"[$([SX1]=[#6])]">(mol);
    validate<"[$([SX1]~P)]">(mol);
    validate<"[$([SX3]=N)]">(mol);
}

template void Rarey_smarts_part_123<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_123<RDKit::ROMol>(RDKit::ROMol &mol);
