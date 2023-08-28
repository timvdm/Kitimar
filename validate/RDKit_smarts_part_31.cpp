#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_31(Mol &mol)
{
    // SMARTS 151 - 155
    validate<"[Si]">(mol);
    validate<"[N;H2]-S">(mol);
    validate<"[N;H]-[C;H]-[N;H]">(mol);
    validate<"c-C(-O)(-O)-c">(mol);
    validate<"c-C(-[O;H])[!O]">(mol);
}

template void RDKit_smarts_part_31<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_31<RDKit::ROMol>(RDKit::ROMol &mol);
