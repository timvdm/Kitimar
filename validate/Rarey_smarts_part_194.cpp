#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_194(Mol &mol)
{
    // SMARTS 966 - 970
    validate<"aC=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"c([OH])c([OH])c([OH])">(mol);
    validate<"c([OH])c([OH])cc([OH])">(mol);
    validate<"c-!@C">(mol);
    validate<"c-!@N">(mol);
}

template void Rarey_smarts_part_194<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_194<RDKit::ROMol>(RDKit::ROMol &mol);
