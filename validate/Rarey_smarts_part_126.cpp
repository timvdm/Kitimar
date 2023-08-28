#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_126(Mol &mol)
{
    // SMARTS 626 - 630
    validate<"[*!H0,#1]">(mol);
    validate<"[+,++,+++]">(mol);
    validate<"[+1]~*~*~[-1]">(mol);
    validate<"[+;!$([+]~[-])]">(mol);
    validate<"[+H]">(mol);
}

template void Rarey_smarts_part_126<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_126<RDKit::ROMol>(RDKit::ROMol &mol);
