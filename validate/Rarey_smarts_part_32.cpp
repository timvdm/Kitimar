#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_32(Mol &mol)
{
    // SMARTS 156 - 160
    validate<"O=C1CCCC(N1)=O">(mol);
    validate<"O=C1NCC2CCCCC21">(mol);
    validate<"O=C1[#6]~[#6]C(=O)[#6]~[#6]1">(mol);
    validate<"O=CN=[N+]=[N-]">(mol);
    validate<"O=[C,N]">(mol);
}

template void Rarey_smarts_part_32<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_32<RDKit::ROMol>(RDKit::ROMol &mol);
