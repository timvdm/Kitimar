#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_51(Mol &mol)
{
    // SMARTS 251 - 255
    validate<"C#C-C(=O)-O-C-C">(mol);
    validate<"C#C-[C;!H0](=O)">(mol);
    validate<"C#C-C(=O)[N,!H0]">(mol);
    validate<"C=C-[C;H1](=O)">(mol);
    validate<"[O;D2]-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
}

template void RDKit_smarts_part_51<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_51<RDKit::ROMol>(RDKit::ROMol &mol);
