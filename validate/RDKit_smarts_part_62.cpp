#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_62(Mol &mol)
{
    // SMARTS 306 - 310
    validate<"[NH;R0][NH;R0]">(mol);
    validate<"[$(O=C[CH](C=O)C=O),$(N#C[CH](-C=O)-C=O)]">(mol);
    validate<"P(=O)(O[H,C])O[H,C]">(mol);
    validate<"[$(N#C-C=[CH][C,c]),$([CH](=[C;R0]-[CH]=O)),$([CH](=[C,R]-C(=O)-C));!$([CH]1=CC(=O)C=CC1=*);!$([CH]1=CC(=O)C(=[N,O])C=C1);!$([CH](=C-C=O)-C=O)]">(mol);
    validate<"[$(N#C-C#C[C,c]),$(C#C-[CH]=O),$(C(#C-C(=O)-[C,c]))]">(mol);
}

template void RDKit_smarts_part_62<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_62<RDKit::ROMol>(RDKit::ROMol &mol);
