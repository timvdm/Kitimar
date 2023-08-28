#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_190(Mol &mol)
{
    // SMARTS 946 - 950
    validate<"[nH1]1ccc(=O)o1">(mol);
    validate<"[nH1]1cnc(n1)C(F)(F)F">(mol);
    validate<"[nH1]1cnnc1C(F)(F)F">(mol);
    validate<"[nH1]1occc1=O">(mol);
    validate<"[nH1]nnn">(mol);
}

template void Rarey_smarts_part_190<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_190<RDKit::ROMol>(RDKit::ROMol &mol);
