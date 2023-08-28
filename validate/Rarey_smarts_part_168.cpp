#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_168(Mol &mol)
{
    // SMARTS 836 - 840
    validate<"[O,N;!H0]-*~*-*=[$([C,N;R0]=O)]">(mol);
    validate<"[O,N]!-CA#A">(mol);
    validate<"[O,N]!-CA=A">(mol);
    validate<"[O;R1][C;R1][C;R1][O;R1][C;R1][C;R1][O;R1]">(mol);
    validate<"[O;X2,X1-][O;X2,X1-]">(mol);
}

template void Rarey_smarts_part_168<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_168<RDKit::ROMol>(RDKit::ROMol &mol);
