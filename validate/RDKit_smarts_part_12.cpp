#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_12(Mol &mol)
{
    // SMARTS 56 - 60
    validate<"[O,N,S][CH2]N1C(=O)CCC1(=O)">(mol);
    validate<"O=C-N!@N">(mol);
    validate<"C!@N=*">(mol);
    validate<"[S,P](=O)(=O)OC">(mol);
    validate<"[$(N!@N),$([N;R0]=[N;R0])]">(mol);
}

template void RDKit_smarts_part_12<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_12<RDKit::ROMol>(RDKit::ROMol &mol);
