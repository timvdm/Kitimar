#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_7(Mol &mol)
{
    // SMARTS 31 - 35
    validate<"C#!@N">(mol);
    validate<"C#@C">(mol);
    validate<"C#C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"C(=*)(=*)">(mol);
    validate<"C(=O)(~c)~c">(mol);
}

template void Rarey_smarts_part_7<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_7<RDKit::ROMol>(RDKit::ROMol &mol);
