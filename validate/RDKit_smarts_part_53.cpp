#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_53(Mol &mol)
{
    // SMARTS 261 - 265
    validate<"C(=O)-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"C(=O)-[C;H2]([F,Br,I,Cl])">(mol);
    validate<"S(=O)(=O)-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"S(=O)(=O)-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"S(=O)(=O)-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
}

template void RDKit_smarts_part_53<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_53<RDKit::ROMol>(RDKit::ROMol &mol);
