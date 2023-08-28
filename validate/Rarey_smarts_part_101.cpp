#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_101(Mol &mol)
{
    // SMARTS 501 - 505
    validate<"[#7]:[#7$(a1aaaa1)]">(mol);
    validate<"[#7]:[#7$(a1aaaaa1)]">(mol);
    validate<"[#7]:[#8$(a1aaaa1)]">(mol);
    validate<"[#7]:[#8$(a1aaaaa1)]">(mol);
    validate<"[#7]=[#6]-1-[#16]-[#6](=[#7])-[#7]=[#6]-1">(mol);
}

template void Rarey_smarts_part_101<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_101<RDKit::ROMol>(RDKit::ROMol &mol);
