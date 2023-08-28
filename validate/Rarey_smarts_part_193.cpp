#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_193(Mol &mol)
{
    // SMARTS 961 - 965
    validate<"a-[I]">(mol);
    validate<"a1aa1">(mol);
    validate<"a1aaaa1">(mol);
    validate<"a1aaaa2aaaa(a12)a1aaaa2aaaaa12">(mol);
    validate<"a1aaaaa1">(mol);
}

template void Rarey_smarts_part_193<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_193<RDKit::ROMol>(RDKit::ROMol &mol);
