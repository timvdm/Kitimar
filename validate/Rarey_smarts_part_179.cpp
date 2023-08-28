#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_179(Mol &mol)
{
    // SMARTS 891 - 895
    validate<"[R](-*(-*))~*~[a]">(mol);
    validate<"[R]~*-[CH3]">(mol);
    validate<"[R]~*~*-[CH3]">(mol);
    validate<"[R]~[D3]">(mol);
    validate<"[S,C](=[O,S])[F,Br,Cl,I]">(mol);
}

template void Rarey_smarts_part_179<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_179<RDKit::ROMol>(RDKit::ROMol &mol);
