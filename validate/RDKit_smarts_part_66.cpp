#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_66(Mol &mol)
{
    // SMARTS 326 - 330
    validate<"C#C">(mol);
    validate<"[Cl,Br,I][CH]C=C">(mol);
    validate<"cCC(=O)[OH]">(mol);
    validate<"C([Cl,Br,I])([Cl,Br,I])([Cl,Br,I])C(=O)[OH]">(mol);
    validate<"[#6]C(=O)C(=O)[OH]">(mol);
}

template void RDKit_smarts_part_66<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_66<RDKit::ROMol>(RDKit::ROMol &mol);
