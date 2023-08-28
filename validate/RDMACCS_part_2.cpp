#include "Validate.hpp"

template<Molecule::Molecule Mol>
void RDMACCS_part_2(Mol &mol)
{
    // SMARTS 6 - 10
    validate<"[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]">(mol);
    validate<"[!#6;!#1]~1~*~*~*~1">(mol);
    validate<"[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]">(mol);
    validate<"[Be,Mg,Ca,Sr,Ba,Ra]">(mol);
    validate<"*~1~*~*~*~1">(mol);
}

template void RDMACCS_part_2<OpenBabel::OBMol>(OpenBabel::OBMol &mol);
template void RDMACCS_part_2<RDKit::ROMol>(RDKit::ROMol &mol);
