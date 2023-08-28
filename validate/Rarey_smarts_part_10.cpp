#include "Validate.hpp"

template<Molecule::Molecule Mol>
void Rarey_smarts_part_10(Mol &mol)
{
    // SMARTS 46 - 50
    validate<"C(=O)[SH1]">(mol);
    validate<"C(=[O,S])(N)Oc">(mol);
    validate<"C(C)(C)=CC=[O,SX2]">(mol);
    validate<"C([F,Cl,Br,I])[CH2,CH][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol);
    validate<"C-!@C">(mol);
}

template void Rarey_smarts_part_10<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void Rarey_smarts_part_10<RDKit::ROMol>(RDKit::ROMol &mol);
