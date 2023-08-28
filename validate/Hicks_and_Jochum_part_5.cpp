#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Hicks_and_Jochum_part_5(Mol &mol)
{
    // SMARTS 21 - 25
    validate<"[#6]=,:1[#6]=,:[#6][#8D2][#6]=,:1">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]=,:1[#6D2]=,:[#6D2][#7D2][#6D3]=,:1-[*;D1,D2,D3,D4]">(mol);
    validate<"[#6]=,:1[#6]=,:[#6][#7D2][#6]=,:1">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]=,:1[#6D2]=,:[#6D2][#16D2][#6D3]=,:1-[*;D1,D2,D3,D4]">(mol);
    validate<"[#6]=,:1[#6]=,:[#6][#16D2][#6]=,:1">(mol);
}

template void Hicks_and_Jochum_part_5<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Hicks_and_Jochum_part_5<RDKit::ROMol>(RDKit::ROMol &mol);
