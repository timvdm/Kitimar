#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_18(Mol &mol)
{
    // SMARTS 86 - 90
    validate<"C=@N">(mol);
    validate<"C=C([F,Cl,Br,I])[F,Cl,Br,I]">(mol);
    validate<"C=CC=CC=CC=C">(mol);
    validate<"C=C[Cl]">(mol);
    validate<"C=C[F]">(mol);
}

template void Rarey_smarts_part_18<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_18<RDKit::ROMol>(RDKit::ROMol &mol);
