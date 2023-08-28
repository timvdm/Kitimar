#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Hicks_and_Jochum_part_4(Mol &mol)
{
    // SMARTS 16 - 20
    validate<"[#6D1]-[#7D3]-1-[#6D3]2=,:[#6D2][#6D2]=,:[#6D2][#6D2]=,:[#6D3]2-[#7]-[#6D3]=,:3[#6D2]=,:[#6D2][#6D2]=,:[#6D2][#6D3]-1=,:3">(mol);
    validate<"[#6D1]-[#7D3]-1-[#6]2=,:[#6][#6]=,:[#6][#6]=,:[#6]2-[#7]-[#6]=,:3[#6]=,:[#6][#6]=,:[#6][#6]-1=,:3">(mol);
    validate<"[#6D1]=[#6D2]-[#6D3]-1=[#6D3](-[#6D1])-[#6D3]=,:2[#6D2]=,:[#6D3]3[#7D2][#6D3](=,:[#6D2][#6D3]-4=,:[#7D2][#6D3](-[#6]=[#6D3]-4-[#6D1])=,:[#6D2][#6D3]5=,:[#6][#6D3](-[#6D1])=,:[#6D3]([#6D2]=,:[#6D3]-1[#7D2]=,:2)[#7D2]5)[#6D3](-[#6D1])=,:[#6D3]3-[#6D2]=[#6D1]">(mol);
    validate<"[#6]-1=[#6]-[#6D3]2=,:[#6][#6D3]3=,:[#6][#6]=,:[#6D3]([#6]=,:[#6D3]4-[#6]=[#6]-[#6D3]([#6]=,:[#6D3]5[#6]=,:[#6][#6D3](=,:[#6][#6D3]-1=,:[#7D2]2)[#7D2]5)=,:[#7D2]4)[#7D2]3">(mol);
    validate<"[*;D1,D2,D3,D4]-[#6D3]=,:1[#6D2]=,:[#6D2][#8D2][#6D3]=,:1-[*;D1,D2,D3,D4]">(mol);
}

template void Hicks_and_Jochum_part_4<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Hicks_and_Jochum_part_4<RDKit::ROMol>(RDKit::ROMol &mol);
