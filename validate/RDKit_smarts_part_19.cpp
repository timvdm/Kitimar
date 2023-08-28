#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_19(Mol &mol)
{
    // SMARTS 91 - 95
    validate<"C1-C-C(=O)-O1">(mol);
    validate<"C=C(-S)-S(=O)">(mol);
    validate<"[C;!R]=[C;!R]-[C;!R](-O)">(mol);
    validate<"C=C-C(=N)">(mol);
    validate<"[C;!R]=[C;!R]-[C;!R](=O)">(mol);
}

template void RDKit_smarts_part_19<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_19<RDKit::ROMol>(RDKit::ROMol &mol);
