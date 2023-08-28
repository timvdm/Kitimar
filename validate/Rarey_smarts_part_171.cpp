#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_171(Mol &mol)
{
    // SMARTS 851 - 855
    validate<"[OH1]C1=NC(=O)CC1=O">(mol);
    validate<"[OH1]C1=NC(=O)NO1">(mol);
    validate<"[OH1]C1=NC(=O)ON1">(mol);
    validate<"[OH1]C1=NC(=O)c2ccccc21">(mol);
    validate<"[OH1]C1=NC=NO1">(mol);
}

template void Rarey_smarts_part_171<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_171<RDKit::ROMol>(RDKit::ROMol &mol);
