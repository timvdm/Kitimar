#include "Validate.hpp"

void Rarey_smarts_part_5(OpenBabel::OBMol &mol)
{
    // SMARTS 201 - 250
    validate_contains<"[!#1]">(mol);
    validate_contains<"[!#1]:1:[!#1]-2:[!#1](:[!#1]:[!#1]:[!#1]:1)-[#7](-[#1])-[#7](-[#6]-2=[#8])-[#6]">(mol);
    validate_contains<"[!#1]:1:[!#1]:[!#1]:[!#1](:[!#1]:[!#1]:1)-[#6](-[#1])=[#6](-[#1])-[#6](-[#7](-[#1])-[#7](-[#1])-c2nnnn2-[#6])=[#8]">(mol);
    validate_contains<"[!#1]:1:[!#1]:[!#1]:[!#1](:[!#1]:[!#1]:1)-[#6](-[#1])=[#6](-[#1])-[#6](-[#7]-c:2:c:c:c:3:c(:c:2):c:c:c(:n:3)-[#7](-[#6])-[#6])=[#8]">(mol);
    //validate_contains<"[!#1]:[!#1]-[#6](-[$([#1]),$([#6]#[#7])])=[#6]-1-[#6]=,:[#6]-[#6](=[$([#8]),$([#7;!R])])-[#6]=,:[#6]-1">(mol); // FIXME: recursive
    validate_contains<"[!#1]:[#6]-[#6](=[#16])-[#7](-[#1])-[#7](-[#1])-[#6]:[!#1]">(mol);
    validate_contains<"[!#1]:[#6]-[#6]-1=[#6](-[#1])-[#6](=[#6](-[#6]#[#7])-[#6](=[#8])-[#7]-1-[#1])-[#6]:[#8]">(mol);
    validate_contains<"[!#6&!#1]=[#6]-1-[#6]=,:[#6]-[#6](=[!#6&!#1])-[#6]=,:[#6]-1">(mol);
    validate_contains<"[!#6]-[CH3]">(mol);
    validate_contains<"[!#6]~*~*~*~[D3]">(mol);
    validate_contains<"[!#6]~*~*~[R]">(mol);
    validate_contains<"[!#6]~*~[R]">(mol);
    //validate_contains<"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]">(mol); // FIXME: recursive
    //validate_contains<"[!$(C=C);!$(C#C)]C(=[O,S,N])[O,S][$(C=C),$(C=N),$(C#C),$(C#N)]">(mol); // FIXME: recursive
    //validate_contains<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[F,Cl,Br,I]">(mol); // FIXME: recursive
    //validate_contains<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N]C(=[O,SX1,N])">(mol); // FIXME: recursive
    //validate_contains<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N][CX4,O,S][$(C=O),a,$(C=C),$(C#C),$(C=N),$(C#N)]">(mol); // FIXME: recursive
    //validate_contains<"[!$(C=C);!$(C#C)]C(=[O,SX1,N])[O,S,N][a]">(mol); // FIXME: recursive
    //validate_contains<"[!$([#6+0]);!$(C(F)(F)F);!$(c(:[!c]):[!c])!$([#6]=,#[!#6])]">(mol); // FIXME: recursive
    //validate_contains<"[!$([#6,F,Cl,Br,I,o,s,nX3,#7v5,#15v5,#16v4,#16v6,*+1,*+2,*+3])]">(mol); // FIXME: recursive
    //validate_contains<"[!$([#6,H0,-,-2,-3])]">(mol); // FIXME: recursive
    validate_contains<"[!H0;#7,#8,#9]">(mol);
    //validate_contains<"[!H0;F,Cl,Br,I,N+,$([OH]-*=[!#6]),+]">(mol); // FIXME: recursive
    validate_contains<"[!H1]">(mol);
    validate_contains<"[!H]">(mol);
    validate_contains<"[!R]">(mol);
    validate_contains<"[!a](=*)~*~*~*~[R]">(mol);
    //validate_contains<"[#15]:[#15$(a1aaaa1)]">(mol); // FIXME: recursive
    //validate_contains<"[#15]:[#15$(a1aaaaa1)]">(mol); // FIXME: recursive
    validate_contains<"[#16!H0]">(mol);
    //validate_contains<"[#16;X2]-1-[#6]=[#6](-[#6]#[#7])-[#6](-[#6])(-[#6]=[#8])-[#6](=[#6]-1-[#7](-[#1])-[#1])-[$([#6]=[#8]),$([#6]#[#7])]">(mol); // FIXME: recursive
    validate_contains<"[#16X2H0]">(mol);
    validate_contains<"[#16X2H0][!#16]">(mol);
    validate_contains<"[#16X2H0][#16X2H0]">(mol);
    validate_contains<"[#16X2H]">(mol);
    validate_contains<"[#16X2][OX2H,OX1H0-]">(mol);
    validate_contains<"[#16X2][OX2H0]">(mol);
    validate_contains<"[#16](=[#8])(=[#8])(-c:1:c:n(-[#6](-[#1])-[#1]):c:n:1)-[#7](-[#1])-c:2:c:n(:n:c:2)-[#6](-[#1])(-[#1])-[#6]:[#6]-[#8]-[#6](-[#1])-[#1]">(mol);
    validate_contains<"[#16](=[#8])(=[#8])-[#7](-[#1])-c:1:c(:c(:c(:s:1)-[#6]-[#1])-[#6]-[#1])-[#6](=[#8])-[#7]-[#1]">(mol);
    //validate_contains<"[#16]-1-[#6](=!@[#7]-[$([#1]),$([#7](-[#1])-[#6]:[#6])])-[#7](-[$([#1]),$([#6]:[#7]:[#6]:[#6]:[#16])])-[#6](=[#8])-[#6]-1=[#6](-[#1])-[#6]:[#6]-[$([#17]),$([#8]-[#6]-[#1])]">(mol); // FIXME: recursive
    //validate_contains<"[#16]-1-[#6](=[#7]-[#6]:[#6])-[#7](-[$([#1]),$([#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#8]),$([#6]:[#6])])-[#6](=[#8])-[#6]-1=[#6](-[#1])-[$([#6]:[#6]:[#6]-[#17]),$([#6]:[!#6&!#1])]">(mol); // FIXME: recursive
    validate_contains<"[#16]-1-[#6](=[#7]-[#7]-[#1])-[#16]-[#6](=[#7]-[#6]:[#6])-[#6]-1=[#7]-[#6]:[#6]">(mol);
    validate_contains<"[#16]-1-[#6](=[#8])-[#7]-[#6](=[#16])-[#6]-1=[#6](-[#1])-[#6]:[#6]">(mol);
    //validate_contains<"[#16]-1-[#6](=[#8])-[#7]-[#6](=[#8])-[#6]-1=[#6](-[#1])-[$([#6]-[#35]),$([#6]:[#6](-[#1]):[#6](-[F,Cl,Br,I]):[#6]:[#6]-[F,Cl,Br,I]),$([#6]:[#6](-[#1]):[#6](-[#1]):[#6]-[#16]-[#6](-[#1])-[#1]),$([#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-[#8]-[#6](-[#1])-[#1]),$([#6]:1:[#6](-[#6](-[#1])-[#1]):[#7](-[#6](-[#1])-[#1]):[#6](-[#6](-[#1])-[#1]):[#6]:1)]">(mol); // FIXME: recursive
    //validate_contains<"[#16]:[#15$(a1aaaa1)]">(mol); // FIXME: recursive
    //validate_contains<"[#16]:[#15$(a1aaaaa1)]">(mol); // FIXME: recursive
    //validate_contains<"[#16]:[#16$(a1aaaa1)]">(mol); // FIXME: recursive
    //validate_contains<"[#16]:[#16$(a1aaaaa1)]">(mol); // FIXME: recursive
    validate_contains<"[#16]=[#6]-1-[#6]=,:[#6]-[!#6&!#1]-[#6]=,:[#6]-1">(mol);
    validate_contains<"[#16]=[#6]-1-[#7](-[#1])-[#6]=[#6]-[#6]-2=[#6]-1-[#6](=[#8])-[#8]-[#6]-2=[#6]-[#1]">(mol);
}
