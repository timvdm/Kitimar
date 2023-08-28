#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_23(Mol &mol)
{
    // SMARTS 111 - 115
    validate<"N#CC(=O)">(mol);
    validate<"N#CC[OH]">(mol);
    validate<"N-!@N">(mol);
    validate<"N-!@O">(mol);
    validate<"N-!@P">(mol);
}

template void Rarey_smarts_part_23<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_23<RDKit::ROMol>(RDKit::ROMol &mol);
