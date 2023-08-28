#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_40(Mol &mol)
{
    // SMARTS 196 - 200
    validate<"N-C(=O)-C(=O)-N">(mol);
    validate<"c=N">(mol);
    validate<"C1-O-C-O-C1">(mol);
    validate<"C=N-N-S(=O)(=O)">(mol);
    validate<"[C;!R](=O)-[C;!R](=O)">(mol);
}

template void RDKit_smarts_part_40<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_40<RDKit::ROMol>(RDKit::ROMol &mol);
