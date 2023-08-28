#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_125(Mol &mol)
{
    // SMARTS 621 - 625
    validate<"[$([nX2r5]:[a-]),$([nX2r5]:[a]:[a-])]">(mol);
    validate<"[$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)]">(mol);
    validate<"[$([nr5]:[nr5,or5,sr5]),$([nr5]:[cr5]:[nr5,or5,sr5])]">(mol);
    validate<"[$(c1(=O)ccn([C,c])cc1),$(c1(=O)n([C,c])cccc1)]">(mol);
    validate<"[$(c:cCl),$(c:c:cCl),$(c:c:c:cCl)]-[$(c:cCl),$(c:c:cCl),$(c:c:c:cCl)]">(mol);
}

template void Rarey_smarts_part_125<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_125<RDKit::ROMol>(RDKit::ROMol &mol);
