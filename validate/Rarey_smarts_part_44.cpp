#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_44(Mol &mol)
{
    // SMARTS 216 - 220
    validate<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N]C(=[O,SX1,N])">(mol);
    validate<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N][CX4,O,S][$(C=O),a,$(C=C),$(C#C),$(C=N),$(C#N)]">(mol);
    validate<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N][a]">(mol);
    validate<"[!$([#6+0]);!$(C(F)(F)F);!$(c(:[!c]):[!c])!$([#6]=,#[!#6])]">(mol);
    validate<"[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]">(mol);
}

template void Rarey_smarts_part_44<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_44<RDKit::ROMol>(RDKit::ROMol &mol);
