#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_173(Mol &mol)
{
    // SMARTS 861 - 865
    validate<"[OH1]C1NC(=O)C(=O)C1">(mol);
    validate<"[OH1]C=C">(mol);
    validate<"[OH1]NC(=O)">(mol);
    validate<"[OH1]NC=[O,S]">(mol);
    validate<"[OH1][C,c,N;!$(C=O)]">(mol);
}

template void Rarey_smarts_part_173<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_173<RDKit::ROMol>(RDKit::ROMol &mol);
