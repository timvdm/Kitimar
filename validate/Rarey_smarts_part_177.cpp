#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_177(Mol &mol)
{
    // SMARTS 881 - 885
    validate<"[OX3H2+]">(mol);
    validate<"[P,S]-;!@[O;D2;$(O([P,S])[#6;!D1][!D1]);!$(OCO);!$(OC=*)]">(mol);
    validate<"[P,S]-;!@[O;D2;$(O[C,c]);!$(OCO);!$(OC=*)]">(mol);
    validate<"[P,S][Cl,Br,F,I]">(mol);
    validate<"[R0;D2][R0;D2][R0;D2][R0;D2]">(mol);
}

template void Rarey_smarts_part_177<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_177<RDKit::ROMol>(RDKit::ROMol &mol);
