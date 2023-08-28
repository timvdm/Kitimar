#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_46(Mol &mol)
{
    // SMARTS 226 - 230
    validate<"[!R]">(mol);
    validate<"[!a](=*)~*~*~*~[R]">(mol);
    validate<"[#15]:[#15$(a1aaaa1)]">(mol);
    validate<"[#15]:[#15$(a1aaaaa1)]">(mol);
    validate<"[#16!H0]">(mol);
}

template void Rarey_smarts_part_46<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_46<RDKit::ROMol>(RDKit::ROMol &mol);
