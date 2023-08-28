#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_8(Mol &mol)
{
    // SMARTS 36 - 40
    validate<"C(=O)-N-C=O">(mol);
    validate<"C(=O)C(=O)">(mol);
    validate<"C(=O)C[N+,n+]">(mol);
    validate<"C(=O)N(C(=O))OC(=O)">(mol);
    validate<"C(=O)OC(=O)">(mol);
}

template void Rarey_smarts_part_8<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_8<RDKit::ROMol>(RDKit::ROMol &mol);
