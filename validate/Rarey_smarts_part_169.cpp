#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_169(Mol &mol)
{
    // SMARTS 841 - 845
    validate<"[OD1H0]=[CD3H0,r6]-[CD3H1,r6]=[CD3H1,r6]-[CD3H0,r6]=[OD1H0]">(mol);
    validate<"[OD1H0]=[CD3H0,r6]-[CD3H1,r6]=[OD1H0]">(mol);
    validate<"[OD1H1]aa[OD1H1]">(mol);
    validate<"[OD1H1]aaaa[OD1H1]">(mol);
    validate<"[OD2H0]-[CD3H0](=[OD1H0])(-[Cl,Br,F])">(mol);
}

template void Rarey_smarts_part_169<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_169<RDKit::ROMol>(RDKit::ROMol &mol);
