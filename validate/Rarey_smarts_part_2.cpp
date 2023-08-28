#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_2(Mol &mol)
{
    // SMARTS 6 - 10
    validate<"*-!:[a]:*:[a]-!:*">(mol);
    validate<"*-!:aa-!:*">(mol);
    validate<"*-!:aaa-!:*">(mol);
    validate<"*-!:aaaa-!:*">(mol);
    validate<"*-!@*">(mol);
}

template void Rarey_smarts_part_2<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_2<RDKit::ROMol>(RDKit::ROMol &mol);
