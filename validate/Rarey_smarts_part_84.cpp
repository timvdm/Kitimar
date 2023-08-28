#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_84(Mol &mol)
{
    // SMARTS 416 - 420
    validate<"[#6]OP(=O)O[#6]">(mol);
    validate<"[#6]S(=O)(=O)O[#6]">(mol);
    validate<"[#6][$([NX2]=O),$(N=C=O),$(OC#N),$(SC#N),$(N=C=S)]">(mol);
    validate<"[#6][CX3](=O)[#6]">(mol);
    validate<"[#6][CX3](=O)[OX2H0][#6]">(mol);
}

template void Rarey_smarts_part_84<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_84<RDKit::ROMol>(RDKit::ROMol &mol);
