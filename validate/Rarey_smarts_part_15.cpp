#include "Validate.hpp"

void Rarey_smarts_part_15(OpenBabel::OBMol &mol)
{
    // SMARTS 701 - 750
    validate_contains<"[CD2H1]=[ND2H0][ND2H1]">(mol);
    validate_contains<"[CD2H1]=[ND2H0][OX2]">(mol);
    validate_contains<"[CD2H1]=[OD1H0]">(mol);
    validate_contains<"[CD2H2]-[CD3H0]=[ND1H1]">(mol);
    validate_contains<"[CD3H0,R](=[SD1H0])([ND2H1,R])([ND2H1,R])">(mol);
    validate_contains<"[CD3H0](-[NX3])(-[NX3])=[SD1H0]">(mol);
    validate_contains<"[CD3H0](-[NX3])(=[NX2])-[SD1H1]">(mol);
    validate_contains<"[CD3H0](=[OD1H0])(-[Cl,Br,F])">(mol);
    validate_contains<"[CD3H0](=[OD1H0])-[OD2H0]-[CD3H0]=[OD1H0]">(mol);
    validate_contains<"[CD3H0]1(=[OD1H0])[NX3][CD3H0](=[OD1H0])[NX3][CD3H0](=[OD1H0])[CX4]1">(mol);
    validate_contains<"[CD3H0]=[CD2H0]=[OD1H0]">(mol);
    validate_contains<"[CD3H1;r6;R2,R3]12-aa-[CD3H1;r6;R2,R3](-aa1)aa2">(mol);
    validate_contains<"[CD3H1]([CD3H0](~[OD1])(~[OD1]))([ND1H2])[CD2H2](a1aaaa2aaaaa12)">(mol);
    //validate_contains<"[CH,CH2,CH3;!$([CH2]CC=[O,S])][F,Cl,Br,I,$(OS(=O)(=O)[#6,#1]),$(OS(=O)(=O)O[#6,#1])]">(mol); // FIXME: recursive
    validate_contains<"[CH2,CH3][NX3][NX2]=[O,S]">(mol);
    validate_contains<"[CH2,CH]=C1C(=[O,SX2])**C1">(mol);
    //validate_contains<"[CH2,CH]=[CH]C=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    //validate_contains<"[CH2,CH]=[CH][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    //validate_contains<"[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1">(mol); // FIXME: recursive
    validate_contains<"[CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]">(mol);
    validate_contains<"[CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]">(mol);
    validate_contains<"[CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]">(mol);
    validate_contains<"[CH2X4][CH2X4][SX2][CH3X4]">(mol);
    validate_contains<"[CH2X4][CHX4]([CH3X4])[CH3X4]">(mol);
    validate_contains<"[CH2X4][CX3](=[OX1])[NX3H2]">(mol);
    validate_contains<"[CH2X4][CX3](=[OX1])[OH0-,OH]">(mol);
    validate_contains<"[CH2X4][OX2H]">(mol);
    validate_contains<"[CH2X4][SX2H,SX1H0-]">(mol);
    validate_contains<"[CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1">(mol);
    validate_contains<"[CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1">(mol);
    validate_contains<"[CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12">(mol);
    //validate_contains<"[CH2]=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    validate_contains<"[CH2]=O">(mol);
    validate_contains<"[CH2]N([CH3])[CH3]">(mol);
    validate_contains<"[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]">(mol);
    validate_contains<"[CH2][NH2]">(mol);
    validate_contains<"[CH3X4]">(mol);
    //validate_contains<"[CHR]=[CR][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    validate_contains<"[CHX4]([CH3X4])[CH2X4][CH3X4]">(mol);
    validate_contains<"[CHX4]([CH3X4])[CH3X4]">(mol);
    validate_contains<"[CHX4]([CH3X4])[OX2H]">(mol);
    validate_contains<"[CH]#C">(mol);
    validate_contains<"[CX1-]#[NX2+]">(mol);
    //validate_contains<"[CX3;$([C]([#6])[#6]),$([CH][#6])]=[NX2][#6]">(mol); // FIXME: recursive
    validate_contains<"[CX3H1](=O)[#6]">(mol);
    validate_contains<"[CX3](=O)[O-]">(mol);
    validate_contains<"[CX3](=O)[OX1H0-,OX2H1]">(mol);
    validate_contains<"[CX3](=O)[OX2H1]">(mol);
    validate_contains<"[CX3](=[OX1])(O)O">(mol);
    validate_contains<"[CX3](=[OX1])([OX2])[OX2H,OX1H0-1]">(mol);
}
