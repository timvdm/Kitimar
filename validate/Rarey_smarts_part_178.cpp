#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_178(Mol &mol)
{
    // SMARTS 886 - 890
    validate<"[R0;D2]~[R0;D2]~[R0;D2]~[R0;D2]">(mol);
    validate<"[R0]">(mol);
    validate<"[R]">(mol);
    validate<"[R](-*(-*))~*~*~*~[a]">(mol);
    validate<"[R](-*(-*))~*~*~[a]">(mol);
}

template void Rarey_smarts_part_178<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_178<RDKit::ROMol>(RDKit::ROMol &mol);
