#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_8(Mol &mol)
{
    // SMARTS 36 - 40
    validate<"[#6]-1=[#6]-[#6D3]2=,:[#6][#6D3]3=,:[#6][#6]=,:[#6D3]([#6D2]=,:[#6D3]4-[#6]=[#6]-[#6D3]([#6]=,:[#6D3]5[#6]=,:[#6][#6D3](=,:[#6][#6D3]-1=,:[#7D2]2)[#7D2]5)=,:[#7D2]4)[#7D2]3">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]=,:1[#6D2]=,:[#6D2][#8D2][#6D3]=,:1-[*;D1,D2,D3,D4]">(mol);
    validate<"[#6]=,:1[#6]=,:[#6][#8D2][#6]=,:1">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]=,:1[#6D2]=,:[#6D2][#7D2][#6D3]=,:1-[*;D1,D2,D3,D4]">(mol);
    validate<"[#6]=,:1[#6]=,:[#6][#7D2][#6]=,:1">(mol);
}

template void Agrafiotis_ABCD_part_8<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_8<RDKit::ROMol>(RDKit::ROMol &mol);
