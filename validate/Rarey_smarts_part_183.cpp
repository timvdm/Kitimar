#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_183(Mol &mol)
{
    // SMARTS 911 - 915
    validate<"[a][CH]=O">(mol);
    validate<"[a]~*~*-[CH3]">(mol);
    validate<"[c,n;R][O;!R][C;!R](=O)[#6,#7,#8;!R]">(mol);
    validate<"[c;$([*Cl]),$([*H1])]1ccc(O)c(C)c1">(mol);
    validate<"[c;$(c([cH])([cH,n&H0,O]))][F,Cl,Br,I]">(mol);
}

template void Rarey_smarts_part_183<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_183<RDKit::ROMol>(RDKit::ROMol &mol);
