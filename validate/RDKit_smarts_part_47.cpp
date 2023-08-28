#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_47(Mol &mol)
{
    // SMARTS 231 - 235
    validate<"C=C-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"C=C-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"C=C-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"C=C-[C;H2]([F,Br,I,Cl])">(mol);
    validate<"CN=C=O">(mol);
}

template void RDKit_smarts_part_47<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_47<RDKit::ROMol>(RDKit::ROMol &mol);
