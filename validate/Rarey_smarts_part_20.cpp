#include "Validate.hpp"

void Rarey_smarts_part_20(OpenBabel::OBMol &mol)
{
    // SMARTS 951 - 1000
    validate_contains<"[nH]1cnoc1=O">(mol);
    validate_contains<"[nX2r5]">(mol);
    validate_contains<"[nX3r5+]:c:n">(mol);
    validate_contains<"[oX2r5]">(mol);
    validate_contains<"[r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25]">(mol);
    validate_contains<"[r;!r3;!r4;!r5;!r6;!r7]">(mol);
    validate_contains<"[sX2r5]">(mol);
    validate_contains<"a-[Br]">(mol);
    validate_contains<"a-[Cl]">(mol);
    validate_contains<"a-[F]">(mol);
    validate_contains<"a-[I]">(mol);
    validate_contains<"a1aa1">(mol);
    validate_contains<"a1aaaa1">(mol);
    validate_contains<"a1aaaa2aaaa(a12)a1aaaa2aaaaa12">(mol);
    validate_contains<"a1aaaaa1">(mol);
    //validate_contains<"aC=C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    validate_contains<"c([OH])c([OH])c([OH])">(mol);
    validate_contains<"c([OH])c([OH])cc([OH])">(mol);
    validate_contains<"c-!@C">(mol);
    validate_contains<"c-!@N">(mol);
    validate_contains<"c-!@O">(mol);
    validate_contains<"c-!@S">(mol);
    validate_contains<"c-@C">(mol);
    validate_contains<"c-@N">(mol);
    validate_contains<"c-@O">(mol);
    validate_contains<"c-@P">(mol);
    validate_contains<"c-@S">(mol);
    validate_contains<"c1(O[CH3])ccc([CH2][CH]=[CH2])cc1">(mol);
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])c([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc1">(mol); // FIXME: recursive
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])c([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cccc1([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])">(mol); // FIXME: recursive
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])ncc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])cc1">(mol); // FIXME: recursive
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])ncc([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])nc1">(mol); // FIXME: recursive
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])ncccc1([F,Cl,Br,I,$(N(=O)~O),$(C#N),$(C=O),$(C(F)(F)F),$(S=O)])">(mol); // FIXME: recursive
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])ncccn1">(mol); // FIXME: recursive
    //validate_contains<"c1([F,Cl,Br,I,$(N(=O)~O)])ncncc1">(mol); // FIXME: recursive
    validate_contains<"c1([OH])c(O[CH3])cccc1">(mol);
    validate_contains<"c1([OH])ccc(O[CH3])cc1">(mol);
    validate_contains<"c1([OH])ccc([CH2][CH]=[CH2])cc1">(mol);
    validate_contains<"c1(c(c(c(n1-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#1])-[#1])-[#1])-[#6](-[#1])-[#1])-[#1])-[#6](=[#8])-[#8]-[#1]">(mol);
    validate_contains<"c1(c-2c(c(n1-[#6](-[#8])=[#8])-[#6](-[#1])-[#1])-[#16]-[#6](-[#1])(-[#1])-[#16]-2)-[#6](-[#1])-[#1]">(mol);
    validate_contains<"c1(coc(c1-[#1])-[#6](=[#16])-[#7]-2-[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[!#1]-[#6](-[#1])(-[#1])-[#6]-2(-[#1])-[#1])-[#1]">(mol);
    //validate_contains<"c1(nn(c(c1-[$([#1]),$([#6]-[#1])])-[#8]-[#1])-c:2:c(:c(:c(:c(:c:2-[#1])-[#1])-[#1])-[#1])-[#1])-[#6;X4]">(mol); // FIXME: recursive
    validate_contains<"c12ccccc1cccc2">(mol);
    validate_contains<"c12ccccc1ccn2">(mol);
    validate_contains<"c1c(-[#7](-[#1])-[#1])nnc1-c2c(-[#6](-[#1])-[#1])oc(c2-[#1])-[#1]">(mol);
    validate_contains<"c1c(=[O,NH2,NH])c(=[O,NH2,NH])ccc1">(mol);
    validate_contains<"c1c(=[O,NH2,NH])ccc(=[O,NH2,NH])c1">(mol);
    //validate_contains<"c1c([OH,NH2,NH])c([OH,NH2,NH,$(N=N),$(N(C)C)])ccc1">(mol); // FIXME: recursive
    //validate_contains<"c1c([OH,NH2,NH])cc([OH,NH2,NH,$(N(C)C)])cc1">(mol); // FIXME: recursive
    //validate_contains<"c1c([OH,NH2,NH])ccc([OH,NH2,NH,$(N=N),$(N(C)C)])c1">(mol); // FIXME: recursive
}
