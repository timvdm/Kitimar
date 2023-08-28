#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_16(Mol &mol)
{
    // SMARTS 76 - 80
    validate<"C1C(=[O,S])[O,S][CH2,CH]1">(mol);
    validate<"C1[O,N]C1">(mol);
    validate<"C1[O,S,N]C1">(mol);
    validate<"C=!@C">(mol);
    validate<"C=!@N">(mol);
}

template void Rarey_smarts_part_16<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_16<RDKit::ROMol>(RDKit::ROMol &mol);
