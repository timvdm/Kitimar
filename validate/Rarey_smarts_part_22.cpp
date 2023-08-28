#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_22(Mol &mol)
{
    // SMARTS 106 - 110
    validate<"C~C~O">(mol);
    validate<"F[CH2,CH]C(F)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"I[CH2,CH]C(I)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"N!@[C$(C@*)]">(mol);
    validate<"N!@[N$(N@*)]">(mol);
}

template void Rarey_smarts_part_22<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_22<RDKit::ROMol>(RDKit::ROMol &mol);
