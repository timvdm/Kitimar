#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_5(Mol &mol)
{
    // SMARTS 21 - 25
    validate<"*1**1">(mol);
    validate<"A-[Cl,Br,I]">(mol);
    validate<"A[#6;!$(C(F)(F)F)]">(mol);
    validate<"Br[CH2,CH]C(Br)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"C!@[C$(C@*)]">(mol);
}

template void Rarey_smarts_part_5<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_5<RDKit::ROMol>(RDKit::ROMol &mol);
