#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_128(Mol &mol)
{
    // SMARTS 636 - 640
    validate<"[C+,Cl+,I+,P+,S+]">(mol);
    validate<"[C,N,O,S,a]c1[n,o,s]c2ccccc2[n,o,s]1">(mol);
    validate<"[C,O]=[#6]1[#7,#8,#16][#6](=[O,N,SX1])c2ccccc12">(mol);
    validate<"[C,P;H1](=[O,S])[O,S]">(mol);
    validate<"[C,c](=N)N">(mol);
}

template void Rarey_smarts_part_128<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_128<RDKit::ROMol>(RDKit::ROMol &mol);
