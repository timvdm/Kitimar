#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_45(Mol &mol)
{
    // SMARTS 221 - 225
    validate<"[!$([#6,H0,-,-2,-3])]">(mol);
    validate<"[!H0;#7,#8,#9]">(mol);
    validate<"[!H0;F,Cl,Br,I,N+,$([OH]-*=[!#6]),+]">(mol);
    validate<"[!H1]">(mol);
    //validate<"[!H]">(mol); // implementation defined
}

template void Rarey_smarts_part_45<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_45<RDKit::ROMol>(RDKit::ROMol &mol);
