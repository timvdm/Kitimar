#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_85(Mol &mol)
{
    // SMARTS 421 - 425
    validate<"[$([#6]);!$(C=[O,S,N])]C(=S)O[$([#6]);!$(C=[O,S,N])]">(mol);
    validate<"[$([#16;D2]([$([CX4,c]);!$(C[O,N,S])])[$([CX4,c]);!$(C[O,N,S])]);!$([#16;D2]1[CX4][CX4]1);!$(s1aaaa1)]">(mol);
    validate<"[CX4][SH]">(mol);
    validate<"c[SH]">(mol);
    validate<"NC(=S)N">(mol);
}

template void RDKit_smarts_part_85<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_85<RDKit::ROMol>(RDKit::ROMol &mol);
