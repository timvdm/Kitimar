#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_144(Mol &mol)
{
    // SMARTS 716 - 720
    validate<"[CH2,CH]=C1C(=[O,SX2])**C1">(mol);
    validate<"[CH2,CH]=[CH]C=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"[CH2,CH]=[CH][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1">(mol);
    validate<"[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]">(mol);
}

template void Rarey_smarts_part_144<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_144<RDKit::ROMol>(RDKit::ROMol &mol);
