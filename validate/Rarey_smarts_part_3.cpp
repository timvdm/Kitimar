#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_3(Mol &mol)
{
    // SMARTS 11 - 15
    validate<"*-*-[R]:[a]">(mol);
    validate<"*-*-[R]~*:[a]">(mol);
    validate<"*-*-[R]~*~*:[a]">(mol);
    //validate<"*/,\[R]=,:;@[R]/,\*">(mol); // FIXME: stereo
    //validate<"*/,\[R]=;@[R]/,\*">(mol); // FIXME: stereo
}

template void Rarey_smarts_part_3<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_3<RDKit::ROMol>(RDKit::ROMol &mol);
