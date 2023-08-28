#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_167(Mol &mol)
{
    // SMARTS 831 - 835
    validate<"[NX3][CX4][NX3]">(mol);
    validate<"[NX3][NX2]=[*]">(mol);
    validate<"[NX3][NX3]">(mol);
    validate<"[NX4]">(mol);
    validate<"[O,N;!H0;R0]">(mol);
}

template void Rarey_smarts_part_167<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_167<RDKit::ROMol>(RDKit::ROMol &mol);
