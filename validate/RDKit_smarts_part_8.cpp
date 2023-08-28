#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_8(Mol &mol)
{
    // SMARTS 36 - 40
    validate<"C1[O,S,N]C1">(mol);
    validate<"c([OH])cc([OH])c([OH])">(mol);
    validate<"c([OH])c([OH])c([OH])">(mol);
    validate<"N=NC(=S)N">(mol);
    validate<"SC#N">(mol);
}

template void RDKit_smarts_part_8<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_8<RDKit::ROMol>(RDKit::ROMol &mol);
