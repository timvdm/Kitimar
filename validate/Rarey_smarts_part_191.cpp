#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_191(Mol &mol)
{
    // SMARTS 951 - 955
    validate<"[nH]1cnoc1=O">(mol);
    validate<"[nX2r5]">(mol);
    validate<"[nX3r5+]:c:n">(mol);
    validate<"[oX2r5]">(mol);
    validate<"[r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25]">(mol);
}

template void Rarey_smarts_part_191<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_191<RDKit::ROMol>(RDKit::ROMol &mol);
