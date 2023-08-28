#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_11(Mol &mol)
{
    // SMARTS 51 - 55
    validate<"O=C1CCCC(N1)=O">(mol);
    validate<"O1CCCCC1C2CCCO2">(mol);
    validate<"[OH]c1cc([OH])cc2OC(C([OH])Cc21)c3cc([OH])c([OH])cc3">(mol);
    validate<"C12OCCC(O1)CC2">(mol);
    validate<"c-N(=O)~O">(mol);
}

template void RDKit_smarts_part_11<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_11<RDKit::ROMol>(RDKit::ROMol &mol);
