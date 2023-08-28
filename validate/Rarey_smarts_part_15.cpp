#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_15(Mol &mol)
{
    // SMARTS 71 - 75
    validate<"C-C[I]">(mol);
    validate<"C1(=O)NS(=O)(=O)[C,c]=,:[C,c]1">(mol);
    validate<"C1(=O)OCC1">(mol);
    validate<"C12OCCC(O1)CC2">(mol);
    validate<"C1=[C,N][$(S(=O)(=O)),$(C=[N,O]),$(S=O)][C,N]=C1c2ccccc2">(mol);
}

template void Rarey_smarts_part_15<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_15<RDKit::ROMol>(RDKit::ROMol &mol);
