#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_21(Mol &mol)
{
    // SMARTS 101 - 105
    validate<"Cl">(mol);
    validate<"[!#6;!#1;!H0]~*~[CH2]~[!#1]">(mol);
    validate<"*@*(@*)@*">(mol);
    validate<"[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]">(mol);
    validate<"[F,Cl,Br,I]~*(~[!#1])~[!#1]">(mol);
}

template void RDMACCS_part_21<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_21<RDKit::ROMol>(RDKit::ROMol &mol);
