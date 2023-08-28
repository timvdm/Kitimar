#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDKit_smarts_part_57(Mol &mol)
{
    // SMARTS 281 - 285
    validate<"S(O[C,c,N,n])(~O)[C,c,N,n]">(mol);
    validate<"[N;R0](~N)~O">(mol);
    validate<"N(~O)(~O)(~O)-*">(mol);
    validate<"[N+]([O-])(=C)-*">(mol);
    validate<"[!$([C,c]-N(=O)~O);$([!O]~[N;R0]=O)]">(mol);
}

template void RDKit_smarts_part_57<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDKit_smarts_part_57<RDKit::ROMol>(RDKit::ROMol &mol);
