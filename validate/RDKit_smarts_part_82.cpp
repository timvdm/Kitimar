#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_82(Mol &mol)
{
    // SMARTS 406 - 410
    validate<"[PX4](=O)(~O)~O">(mol);
    validate<"O=P(~O)(~O)(~O)">(mol);
    validate<"[CX4][NH2]">(mol);
    validate<"[c,C]1(~[O;D1])~*!-*~[c,C](~[O;D1])~*!-*~1">(mol);
    validate<"[$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]">(mol);
}

template void RDKit_smarts_part_82<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_82<RDKit::ROMol>(RDKit::ROMol &mol);
