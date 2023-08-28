#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_16(Mol &mol)
{
    // SMARTS 76 - 80
    validate<"C(=O)-C=C-C(=O)">(mol);
    validate<"C(=O)-N(-O)-C(=O)">(mol);
    validate<"C(=O)-N(-S)-C(=O)">(mol);
    validate<"[C;!R](=O)-[N;!R]-[C;!R](=O)">(mol);
    validate<"C(=O)-N-N(=O)">(mol);
}

template void RDKit_smarts_part_16<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_16<RDKit::ROMol>(RDKit::ROMol &mol);
