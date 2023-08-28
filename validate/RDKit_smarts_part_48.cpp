#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_48(Mol &mol)
{
    // SMARTS 236 - 240
    validate<"CN=C=S">(mol);
    validate<"C(=O)([F,I,Br,Cl])">(mol);
    validate<"[N;!R]=[N;!R]">(mol);
    validate<"[S;!R]-[S;!R]">(mol);
    validate<"C1NC(=O)OC(=O)C1">(mol);
}

template void RDKit_smarts_part_48<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_48<RDKit::ROMol>(RDKit::ROMol &mol);
