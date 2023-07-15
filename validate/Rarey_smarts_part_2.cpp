#include "Validate.hpp"

void Rarey_smarts_part_2(OpenBabel::OBMol &mol)
{
    // SMARTS 51 - 100
    validate_contains<"C-!@N">(mol);
    validate_contains<"C-!@O">(mol);
    validate_contains<"C-!@P">(mol);
    validate_contains<"C-!@S">(mol);
    validate_contains<"C-!@[Br]">(mol);
    validate_contains<"C-!@[Cl]">(mol);
    validate_contains<"C-!@[F]">(mol);
    validate_contains<"C-!@[I]">(mol);
    validate_contains<"C-@A[Br]">(mol);
    validate_contains<"C-@A[Cl]">(mol);
    validate_contains<"C-@A[F]">(mol);
    validate_contains<"C-@A[I]">(mol);
    validate_contains<"C-@C">(mol);
    validate_contains<"C-@N">(mol);
    validate_contains<"C-@O">(mol);
    validate_contains<"C-@P">(mol);
    validate_contains<"C-@S">(mol);
    validate_contains<"C-C[Br]">(mol);
    validate_contains<"C-C[Cl]">(mol);
    validate_contains<"C-C[F]">(mol);
    validate_contains<"C-C[I]">(mol);
    validate_contains<"C1(=O)NS(=O)(=O)[C,c]=,:[C,c]1">(mol);
    validate_contains<"C1(=O)OCC1">(mol);
    validate_contains<"C12OCCC(O1)CC2">(mol);
    //validate_contains<"C1=[C,N][$(S(=O)(=O)),$(C=[N,O]),$(S=O)][C,N]=C1c2ccccc2">(mol); // FIXME: recursive
    validate_contains<"C1C(=[O,S])[O,S][CH2,CH]1">(mol);
    validate_contains<"C1[O,N]C1">(mol);
    validate_contains<"C1[O,S,N]C1">(mol);
    validate_contains<"C=!@C">(mol);
    validate_contains<"C=!@N">(mol);
    validate_contains<"C=!@O">(mol);
    validate_contains<"C=!@P">(mol);
    validate_contains<"C=!@S">(mol);
    validate_contains<"C=@A[I]">(mol);
    validate_contains<"C=@C">(mol);
    validate_contains<"C=@N">(mol);
    validate_contains<"C=C([F,Cl,Br,I])[F,Cl,Br,I]">(mol);
    validate_contains<"C=CC=CC=CC=C">(mol);
    validate_contains<"C=C[Cl]">(mol);
    validate_contains<"C=C[F]">(mol);
    validate_contains<"C=P">(mol);
    validate_contains<"C=[OD2H0,SD2H0]-[SD1H0,OD1H0]">(mol);
    validate_contains<"C=[SD1H0]">(mol);
    validate_contains<"CC(C)=[CH][CH2][OH]">(mol);
    validate_contains<"CC=C(C)[CH2][OH]">(mol);
    validate_contains<"COS(=O)(=O)[C,c]">(mol);
    validate_contains<"COS(=O)O[C,c]">(mol);
    //validate_contains<"C[C@?H](Cl)Br">(mol); // FIXME: stereo
    //validate_contains<"C[C@?](F)(Cl)Br">(mol); // FIXME: stereo
    validate_contains<"C[O,S;R0][C;R0](=S)">(mol);
}
