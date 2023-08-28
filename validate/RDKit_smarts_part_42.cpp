#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_42(Mol &mol)
{
    // SMARTS 206 - 210
    validate<"C1(O)CC(O)CO1">(mol);
    validate<"C1(O)C(O)CCO1">(mol);
    validate<"C1(O)CC(O)CCO1">(mol);
    validate<"N=[N+]=N">(mol);
    validate<"N-C#N">(mol);
}

template void RDKit_smarts_part_42<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_42<RDKit::ROMol>(RDKit::ROMol &mol);
