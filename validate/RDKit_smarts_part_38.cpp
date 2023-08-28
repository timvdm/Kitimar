#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_38(Mol &mol)
{
    // SMARTS 186 - 190
    validate<"C([F,Br,Cl,I])-S">(mol);
    validate<"C(-S)(-S)=C(-S)(-S)">(mol);
    validate<"C(=N)-C=N-O">(mol);
    validate<"C(-C#N)(-C#N)">(mol);
    validate<"C-O-N=O">(mol);
}

template void RDKit_smarts_part_38<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_38<RDKit::ROMol>(RDKit::ROMol &mol);
