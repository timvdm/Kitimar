#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_110(Mol &mol)
{
    // SMARTS 546 - 550
    validate<"[#9,#17,#35,#53]-[R]">(mol);
    validate<"[#9,#17,#35,#53]-[a]">(mol);
    validate<"[#9,#17,#35,#53]~*(~*)~*">(mol);
    validate<"[#9,#17,#35,#53]~[D3]">(mol);
    validate<"[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]">(mol);
}

template void Rarey_smarts_part_110<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_110<RDKit::ROMol>(RDKit::ROMol &mol);
