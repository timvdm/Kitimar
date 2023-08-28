#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_21(Mol &mol)
{
    // SMARTS 101 - 105
    validate<"C[OX2][CX3](=[OX1])[OX2]C">(mol);
    validate<"Cc1:c(O):c:c:[$(cCl),$([cH])]:c1">(mol);
    validate<"Cl[CH2,CH]C(Cl)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"C~C(~O)~O">(mol);
    validate<"C~C~C~C">(mol);
}

template void Rarey_smarts_part_21<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_21<RDKit::ROMol>(RDKit::ROMol &mol);
