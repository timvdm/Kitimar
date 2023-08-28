#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Agrafiotis_ABCD_part_5(Mol &mol)
{
    // SMARTS 21 - 25
    validate<"[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=[#8D1])-[#8]">(mol);
    validate<"[#6D1]-[#6]-1-2-[#6]-[#6D2]-[#6D2]-[#6D3]-1-[#6D3]-3-[#6D2]-[#6D2]-[#6D3]-4-[#6D2]-[#6]-[#6D2]-[#6D2]-[#6]-4(-[#6D1])-[#6D3]-3-[#6D2]-[#6D2]-2">(mol);
    validate<"[#6]-[#6]-1-2-[#6]-[#6]-[#6]-[#6]-1-[#6]-3-[#6]-[#6]-[#6]-4-[#6]-[#6]-[#6]-[#6]-[#6]-4(-[#6])-[#6]-3-[#6]-[#6]-2">(mol);
    validate<"[#6D2]-1=[#6D2]-[#6D2]-2-[#6]-[#6]-[#6D3]-1-[#6D2]-2">(mol);
    validate<"[#6]-1=[#6]-[#6]-2-[#6]-[#6]-[#6]-1-[#6]-2">(mol);
}

template void Agrafiotis_ABCD_part_5<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Agrafiotis_ABCD_part_5<RDKit::ROMol>(RDKit::ROMol &mol);
