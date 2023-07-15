#include "Validate.hpp"

void Rarey_smarts_part_12(OpenBabel::OBMol &mol)
{
    // SMARTS 551 - 600
    //validate_contains<"[$(C(=O)(~c)~c);!$([$(c1(=O)ccn([C,c])cc1),$(c1(=O)n([C,c])cccc1)])]">(mol); // FIXME: recursive
    //validate_contains<"[$(P(=O)[O,S]);!$(P[OH1])]">(mol); // FIXME: recursive
    //validate_contains<"[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]">(mol); // FIXME: recursive
    //validate_contains<"[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]">(mol); // FIXME: recursive
    //validate_contains<"[$([!-0!-1!-2!-3!-4]~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~*~*~[!+0!+1!+2!+3!+4]),$([!-0!-1!-2!-3!-4]~*~*~*~*~*~*~*~*~*~[!+0!+1!+2!+3!+4])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X3](=[OX1])[OX2H,OX1H0-]),$([#16X3+]([OX1-])[OX2H,OX1H0-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X3](=[OX1])[OX2H0]),$([#16X3+]([OX1-])[OX2H0])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X3]=[OX1]),$([#16X3+][OX1-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2H,OX1H0-]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2H,OX1H0-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#16X4]([NX3])(=[OX1])(=[OX1])[OX2][#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[OX2][#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#1X1][$([NX4+]),$([NX3]);!$(*=*)&!$(*:*)])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#1X1][$([nX3](:*):*),$([nX2](:*):*),$([#7X2]=*),$([NX3](=*)=*),$([#7X3+](-*)=*),$([#7X3+H]=*)])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#1]),$([#6](-[#1])-[#1]),$([#6]:[#6])]-c:1:c(:c(:c(:s:1)-[#7](-[#1])-[#6](=[#8])-[#6])-[#6](=[#8])-[#8])-[$([#6]:1:[#6]:[#6]:[#6]:[#6]:[#6]:1),$([#6]:1:[#16]:[#6]:[#6]:[#6]:1)]">(mol); // FIXME: recursive
    //validate_contains<"[$([#1]),$([#6](-[#1])-[#1])]-[#8]-c:1:c(:c(:c(:c(:c:1-[#1])-[#1])-[#1])-[#1])-[#7](-[#1])-[#6](-[#1])(-[#1])-c:2:n:c:c:n:2">(mol); // FIXME: recursive
    //validate_contains<"[$([#6+0]);!$(C(F)(F)F);!$(c(:[!c]):[!c])!$([#6]=,#[!#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#6X4@](*)(*)(*)*),$([#6X4@H](*)(*)*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([#6]-!@[A!C!H]-!@[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]">(mol); // FIXME: recursive
    //validate_contains<"[$([#7](-[#1])-[#1]),$([#8]-[#1])]-[#6]-2=[#6](-[#6]#[#7])-[#6](-[#1])(-[#6]:[#6])-c:1:c(:n(-[#6]):n:c:1)-[#8]-2">(mol); // FIXME: recursive
    //validate_contains<"[$([#8X2v2]([#6])[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N])]">(mol); // FIXME: recursive
    //validate_contains<"[$([$([NX3]=O),$([NX3+][O-])])]">(mol); // FIXME: recursive
    //validate_contains<"[$([$([NX4]=O),$([NX4+][O-])])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CH]=[CH2,CH]),$(C(C)=[CH2,CH]),$(C#C);!$(C(C)=CC)][CH2][OH]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX2]#C)]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX2](=C)=C)]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]=[CX3])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]=[OX1]),$([CX3+]-[OX1-])]">(mol); // FIXME: recursive
    //validate_contains<"[$([CX3]=[OX1])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH2]!:c),$([NH1]([CX4])!:c),$([NH0]([CX4])([CX4])!:c)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NH2]-C),$([OH]-*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX1]#*)]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]);!$([$([$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H2][CX3](=[OX1])[OX2H,OX1-,N]),$([$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])[CX3](=[OX1])[OX2H,OX1-,N])])]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]">(mol); // FIXME: recursive
    //validate_contains<"[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]([$([CH3X4]),$([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3]),$([CH2X4][CX3](=[OX1])[NX3H2]),$([CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][SX2H,SX1H0-]),$([CH2X4][CH2X4][CX3](=[OX1])[OH0-,OH]),$([CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1),$([CHX4]([CH3X4])[CH2X4][CH3X4]),$([CH2X4][CHX4]([CH3X4])[CH3X4]),$([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0]),$([CH2X4][CH2X4][SX2][CH3X4]),$([CH2X4][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1),$([CH2X4][OX2H]),$([CHX4]([CH3X4])[OX2H]),$([CH2X4][cX3]1[cX3H][nX3H][cX3]2[cX3H][cX3H][cX3H][cX3H][cX3]12),$([CH2X4][cX3]1[cX3H][cX3H][cX3]([OHX2,OH0X1-])[cX3H][cX3H]1),$([CHX4]([CH3X4])[CH3X4])])">(mol); // FIXME: recursive
}
