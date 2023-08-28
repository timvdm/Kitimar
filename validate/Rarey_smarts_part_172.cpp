#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_172(Mol &mol)
{
    // SMARTS 856 - 860
    validate<"[OH1]C1=NN=CO1">(mol);
    validate<"[OH1]C1=NOC=C1">(mol);
    validate<"[OH1]C1=NS(=O)(=O)c2ccccc21">(mol);
    validate<"[OH1]C1=NSN=C1">(mol);
    validate<"[OH1]C1=N[NH1]C=N1">(mol);
}

template void Rarey_smarts_part_172<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_172<RDKit::ROMol>(RDKit::ROMol &mol);
