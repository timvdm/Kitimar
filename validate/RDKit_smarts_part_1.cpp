#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_1(Mol &mol)
{
    // SMARTS 1 - 5
    validate<"[Br,Cl,I][CX4;CH,CH2]">(mol);
    validate<"[S,C](=[O,S])[F,Br,Cl,I]">(mol);
    validate<"O=CN=[N+]=[N-]">(mol);
    validate<"COS(=O)O[C,c]">(mol);
    validate<"COS(=O)(=O)[C,c]">(mol);
}

template void RDKit_smarts_part_1<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_1<RDKit::ROMol>(RDKit::ROMol &mol);
