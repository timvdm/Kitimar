#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_52(Mol &mol)
{
    // SMARTS 256 - 260
    validate<"[C;!H3]-O-S(=O)(=O)O">(mol);
    validate<"[C;!a]-S(=O)(=O)[O-]">(mol);
    validate<"[N;!a]-S(=O)(=O)[O-]">(mol);
    validate<"C(=O)-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
    validate<"C(=O)-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
}

template void RDKit_smarts_part_52<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_52<RDKit::ROMol>(RDKit::ROMol &mol);
