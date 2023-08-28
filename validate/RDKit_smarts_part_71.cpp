#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_71(Mol &mol)
{
    // SMARTS 351 - 355
    validate<"cB([OH])[OH]">(mol);
    validate<"[N;R0;D2]~[N;R0]~[N;R0;D1]">(mol);
    validate<"[N;D2]([C,c;!$(C=[O,S,N])])=[N;D2]-[C,c;!$(C=[O,S,N])]">(mol);
    validate<"[OH][CX4][CX4][Cl,Br,I]">(mol);
    validate<"[OH][CX4][CX4][$([NH2]),$([NH][CX4]),$(N([CX4])[CX4])]">(mol);
}

template void RDKit_smarts_part_71<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_71<RDKit::ROMol>(RDKit::ROMol &mol);
