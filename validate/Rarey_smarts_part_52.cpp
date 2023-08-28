#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_52(Mol &mol)
{
    // SMARTS 256 - 260
    validate<"[#6!H0]~@[#6!H0]~@[#6!H0]~@[#6!H0]~@[#6!H0]">(mol);
    validate<"[#6+]">(mol);
    validate<"[#6,#7;R0]=[#8]">(mol);
    //validate<"[#6;!D1;!$([#6]-;!@[#6]=;!@[C,O,S;!@]);$([#6][!H;!D1]!=[!H;!D1])]=;!@[#6;!D1;!$([#6]-;!@[#6]=;!@[C,O,S;!@]);$([#6][!H;!D1]!=[!H;!D1])]">(mol); // FIXME: stereo
    validate<"[#6;X3v3+0]">(mol);
}

template void Rarey_smarts_part_52<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_52<RDKit::ROMol>(RDKit::ROMol &mol);
