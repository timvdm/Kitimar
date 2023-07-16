#include "Validate.hpp"

void Rarey_smarts_part_17(OpenBabel::OBMol &mol)
{
    // SMARTS 801 - 850
    validate_contains<"[NH1]([P,S]=O)([P,S]=O)">(mol);
    validate_contains<"[NH1]1C(=O)c2ccccc2S1(=O)=O">(mol);
    validate_contains<"[NH1]1C=NOS1=O">(mol);
    validate_contains<"[NH1]1[nH]cnc1=O">(mol);
    validate_contains<"[NH2][CX4]">(mol);
    validate_contains<"[NH]([CX4])[CX4]">(mol);
    validate_contains<"[NH](c)S(=O)=O">(mol);
    validate_contains<"[NH]S(=O)(=O)C(F)(F)F">(mol);
    validate_contains<"[NX1]#[CX2]">(mol);
    validate_contains<"[NX2-]">(mol);
    validate_contains<"[NX2]=N">(mol);
    validate_contains<"[NX2]=[NX2]">(mol);
    validate_contains<"[NX2]=[OX1]">(mol);
    validate_contains<"[NX3+]=[CX3]">(mol);
    validate_contains<"[NX3,NX4+][CX3](=[OX1])[OX2,OX1-]">(mol);
    validate_contains<"[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]">(mol);
    validate_contains<"[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]">(mol);
    //validate_contains<"[NX3;H2,H1;!$(NC=O)]">(mol); // FIXME: recursive
    //validate_contains<"[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]">(mol); // FIXME: recursive
    validate_contains<"[NX3H2,NH3X4+][CX4H]([*])[CX3](=[OX1])[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-]">(mol);
    validate_contains<"[NX3]([CX4])([CX4])[CX4]">(mol);
    validate_contains<"[NX3]-[Cl,Br,F,I]">(mol);
    validate_contains<"[NX3]-[OX2]">(mol);
    //validate_contains<"[NX3][$(C=C),$(cc)]">(mol); // FIXME: recursive
    validate_contains<"[NX3][CX2]#[NX1]">(mol);
    validate_contains<"[NX3][CX3](=[OX1])[#6]">(mol);
    validate_contains<"[NX3][CX3](=[OX1])[OX2H0]">(mol);
    validate_contains<"[NX3][CX3]=[CX3]">(mol);
    validate_contains<"[NX3][CX3]=[NX3+]">(mol);
    validate_contains<"[NX3][CX3]=[SX1]">(mol);
    validate_contains<"[NX3][CX4][NX3]">(mol);
    validate_contains<"[NX3][NX2]=[*]">(mol);
    validate_contains<"[NX3][NX3]">(mol);
    validate_contains<"[NX4]">(mol);
    validate_contains<"[O,N;!H0;R0]">(mol);
    //validate_contains<"[O,N;!H0]-*~*-*=[$([C,N;R0]=O)]">(mol); // FIXME: recursive
    validate_contains<"[O,N]!-CA#A">(mol);
    validate_contains<"[O,N]!-CA=A">(mol);
    validate_contains<"[O;R1][C;R1][C;R1][O;R1][C;R1][C;R1][O;R1]">(mol);
    validate_contains<"[O;X2,X1-][O;X2,X1-]">(mol);
    validate_contains<"[OD1H0]=[CD3H0,r6]-[CD3H1,r6]=[CD3H1,r6]-[CD3H0,r6]=[OD1H0]">(mol);
    validate_contains<"[OD1H0]=[CD3H0,r6]-[CD3H1,r6]=[OD1H0]">(mol);
    validate_contains<"[OD1H1]aa[OD1H1]">(mol);
    validate_contains<"[OD1H1]aaaa[OD1H1]">(mol);
    validate_contains<"[OD2H0]-[CD3H0](=[OD1H0])(-[Cl,Br,F])">(mol);
    validate_contains<"[OD2H0][CX4][OD2H0]">(mol);
    validate_contains<"[OD2]([#6])[#6]">(mol);
    validate_contains<"[OH1]C1=CC(=O)NO1">(mol);
    validate_contains<"[OH1]C1=CC(=O)ON1">(mol);
    validate_contains<"[OH1]C1=COC=CC1=O">(mol);
}