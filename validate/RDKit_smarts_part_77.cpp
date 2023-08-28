#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_77(Mol &mol)
{
    // SMARTS 381 - 385
    validate<"[$([NH;R0]([C;!$(C=[O,S,N])]))][$([NH;R0][C;!$(C=[O,S,N])])]">(mol);
    validate<"C=N[NH2]">(mol);
    validate<"C=NC=O">(mol);
    validate<"[$([C;R0]=[N;R0]);!$(C(~[N,O,S])(~[N,O,S]));!$([C;R0]=[N;R0]~[N,O,n])]">(mol);
    validate<"C#N-[#6]">(mol);
}

template void RDKit_smarts_part_77<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_77<RDKit::ROMol>(RDKit::ROMol &mol);
