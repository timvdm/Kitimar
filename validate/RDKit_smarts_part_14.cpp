#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_14(Mol &mol)
{
    // SMARTS 66 - 70
    validate<"O=CC([$(C(F)(F)F),$(C#N),$(Cl)])=C">(mol);
    validate<"C#C-[F,Br,I,Cl]">(mol);
    validate<"C(-[O;H1])(-C#N)">(mol);
    validate<"C(-O)-C-N(=O)=O">(mol);
    validate<"C(=N)-S">(mol);
}

template void RDKit_smarts_part_14<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_14<RDKit::ROMol>(RDKit::ROMol &mol);
