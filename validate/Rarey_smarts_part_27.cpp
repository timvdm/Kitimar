#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_27(Mol &mol)
{
    // SMARTS 131 - 135
    validate<"N=!@O">(mol);
    validate<"N=!@P">(mol);
    validate<"N=@N">(mol);
    validate<"N=C=N">(mol);
    validate<"N=C=[S,O]">(mol);
}

template void Rarey_smarts_part_27<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_27<RDKit::ROMol>(RDKit::ROMol &mol);
