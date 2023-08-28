#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_13(Mol &mol)
{
    // SMARTS 61 - 65
    validate<"[$(N!@S),$(N~C(~S)-N)]">(mol);
    validate<"[CH2]=[CH]-[N,O,S]">(mol);
    validate<"S-C#N">(mol);
    validate<"S=C-[#6,N,O]">(mol);
    validate<"[Cl,Br,I]-N">(mol);
}

template void RDKit_smarts_part_13<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_13<RDKit::ROMol>(RDKit::ROMol &mol);
