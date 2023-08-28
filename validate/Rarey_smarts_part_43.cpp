#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_43(Mol &mol)
{
    // SMARTS 211 - 215
    validate<"[!#6]~*~*~[R]">(mol);
    validate<"[!#6]~*~[R]">(mol);
    validate<"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]">(mol);
    validate<"[!$(C=C);!$(C#C)]C(=[O,S,N])[O,S][$(C=C),$(C=N),$(C#C),$(C#N)]">(mol);
    validate<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[F,Cl,Br,I]">(mol);
}

template void Rarey_smarts_part_43<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_43<RDKit::ROMol>(RDKit::ROMol &mol);
