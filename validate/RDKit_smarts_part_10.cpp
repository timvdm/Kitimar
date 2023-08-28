#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_10(Mol &mol)
{
    // SMARTS 46 - 50
    validate<"N1CCC1=O">(mol);
    validate<"O=C1[#6]~[#6]C(=O)[#6]~[#6]1">(mol);
    validate<"C=CC=CC=CC=C">(mol);
    validate<"O1CCCCC1OC2CCC3CCCCC3C2">(mol);
    validate<"O=C1NCC2CCCCC21">(mol);
}

template void RDKit_smarts_part_10<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_10<RDKit::ROMol>(RDKit::ROMol &mol);
