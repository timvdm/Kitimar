#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_161(Mol &mol)
{
    // SMARTS 801 - 805
    validate<"[NH1]([P,S]=O)([P,S]=O)">(mol);
    validate<"[NH1]1C(=O)c2ccccc2S1(=O)=O">(mol);
    validate<"[NH1]1C=NOS1=O">(mol);
    validate<"[NH1]1[nH]cnc1=O">(mol);
    validate<"[NH2][CX4]">(mol);
}

template void Rarey_smarts_part_161<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_161<RDKit::ROMol>(RDKit::ROMol &mol);
