#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_9(Mol &mol)
{
    // SMARTS 41 - 45
    validate<"cC[N+]">(mol);
    validate<"C[O,S;R0][C;R0](=S)">(mol);
    validate<"N[CH2]C#N">(mol);
    validate<"C1(=O)OCC1">(mol);
    validate<"P(=O)([OH])OP(=O)[OH]">(mol);
}

template void RDKit_smarts_part_9<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_9<RDKit::ROMol>(RDKit::ROMol &mol);
