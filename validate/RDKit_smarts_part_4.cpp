#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_4(Mol &mol)
{
    // SMARTS 16 - 20
    validate<"C(=O)C[N+,n+]">(mol);
    validate<"[N;R0][N;R0]C(=O)">(mol);
    validate<"[C+,Cl+,I+,P+,S+]">(mol);
    validate<"C=P">(mol);
    validate<"[Cl]C([C&R0])=N">(mol);
}

template void RDKit_smarts_part_4<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_4<RDKit::ROMol>(RDKit::ROMol &mol);
