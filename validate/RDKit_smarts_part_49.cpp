#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_49(Mol &mol)
{
    // SMARTS 241 - 245
    validate<"[C;!R]C(=N)S[C;!R]">(mol);
    validate<"[C;!R]C(=N)O[C;!R]">(mol);
    validate<"C=C=C">(mol);
    validate<"C(=O)-O-C(=O)-[!N]">(mol);
    validate<"C=C-C#N">(mol);
}

template void RDKit_smarts_part_49<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_49<RDKit::ROMol>(RDKit::ROMol &mol);
