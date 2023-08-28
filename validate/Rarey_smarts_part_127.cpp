#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_127(Mol &mol)
{
    // SMARTS 631 - 635
    validate<"[+]">(mol);
    validate<"[-,--,---]">(mol);
    validate<"[-;!$([-]~[+])]">(mol);
    validate<"[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]">(mol);
    validate<"[Br,Cl,I][CX4,CH,CH2]">(mol);
}

template void Rarey_smarts_part_127<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_127<RDKit::ROMol>(RDKit::ROMol &mol);
