#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_74(Mol &mol)
{
    // SMARTS 366 - 370
    validate<"ClC(=O)O[CX4,c]">(mol);
    validate<"[CX4,c]-C#N">(mol);
    validate<"N[CX4]C(=O)N[CX4]C(=O)">(mol);
    validate<"[S;D2]-[S;D2]">(mol);
    validate<"[$([S;D2]([CX4,c])!@[CH,CH2]!@[S;D2][CX4,c])]">(mol);
}

template void RDKit_smarts_part_74<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_74<RDKit::ROMol>(RDKit::ROMol &mol);
