#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_1(Mol &mol)
{
    // SMARTS 1 - 5
    validate<"[#103,#104]">(mol);
    validate<"[Ge,#33,#34,Sn,Sb,#52,Tl,Pb,Bi]">(mol);
    validate<"[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]">(mol);
    validate<"[Sc,Ti,Y,Zr,Hf]">(mol);
    validate<"[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]">(mol);
}

template void RDMACCS_part_1<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_1<RDKit::ROMol>(RDKit::ROMol &mol);
