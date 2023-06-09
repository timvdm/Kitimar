#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Util/Util.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>

#include <openbabel/parsmart.h>
#include <gtest/gtest.h>
#include <tuple>

#include "TestData.hpp"


using namespace Kitimar;

template<ctll::fixed_string SMARTS>
void validate_contains(OpenBabel::OBMol &mol)
{
    //std::cout << "    " << Util::toString(SMARTS) <<  std::endl;

    auto ctMatch = CTSmarts::contains<SMARTS>(mol);

    OpenBabel::OBSmartsPattern obsmarts;
    obsmarts.Init(Util::toString(SMARTS));
    auto obMatch = obsmarts.Match(mol);
    if (ctMatch != obMatch)
        std::cout << "FAIL: " << Util::toString(SMARTS) << std::endl;
    EXPECT_EQ(ctMatch, obMatch);
}


TEST(ValidateCTSmarts, SqcRDkit)
{
    OpenBabelSmilesMolSource source{chembl_smi_filename("100K")};
    auto i = 0;
    for (auto mol : source.molecules()) {
        std::cout << "Molecule: " << ++i << " -- " << writeSmiles(mol) << std::endl;

        validate_contains<"[Br,Cl,I][CX4;CH,CH2]">(mol);
        validate_contains<"[S,C](=[O,S])[F,Br,Cl,I]">(mol);
        validate_contains<"O=CN=[N+]=[N-]">(mol);
        validate_contains<"COS(=O)O[C,c]">(mol);
        validate_contains<"COS(=O)(=O)[C,c]">(mol);
        validate_contains<"C(=O)OC(=O)">(mol);
        validate_contains<"OO">(mol);
        validate_contains<"C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)">(mol);
        validate_contains<"C(=O)Oc1ccc(N(=O)=O)cc1">(mol);
        validate_contains<"C(=O)Onnn">(mol);
        validate_contains<"N=C=[S,O]">(mol);
        validate_contains<"OS(=O)(=O)C(F)(F)F">(mol);
        validate_contains<"P(=S)(S)S">(mol);
        validate_contains<"NP(=O)(N)N">(mol);
        validate_contains<"cN=[N+]=[N-]">(mol);
        validate_contains<"C(=O)C[N+,n+]">(mol);
        validate_contains<"[N;R0][N;R0]C(=O)">(mol);
        validate_contains<"[C+,Cl+,I+,P+,S+]">(mol);
        validate_contains<"C=P">(mol);
        validate_contains<"[Cl]C([C&R0])=N">(mol);
        validate_contains<"[N&D2](=O)">(mol);
        validate_contains<"[P,S][Cl,Br,F,I]">(mol);
        validate_contains<"N=C=N">(mol);
        validate_contains<"[N+]#[C-]">(mol);
        validate_contains<"C(=O)N(C(=O))OC(=O)">(mol);
        validate_contains<"N#CC[OH]">(mol);
        validate_contains<"N#CC(=O)">(mol);
        validate_contains<"S(=O)(=O)C#N">(mol);
        validate_contains<"P(OCC)(OCC)(=O)C#N">(mol);
        validate_contains<"[N;R0]=[N;R0]C#N">(mol);
        validate_contains<"[N;R0]=[N;R0]CC=O">(mol);
        validate_contains<"[CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0][CD2;R0]">(mol);
        validate_contains<"[O;R1][C;R1][C;R1][O;R1][C;R1][C;R1][O;R1]">(mol);
        validate_contains<"SS">(mol);
        validate_contains<"[SH]">(mol);
        validate_contains<"C1[O,S,N]C1">(mol);
        validate_contains<"c([OH])cc([OH])c([OH])">(mol);
        validate_contains<"c([OH])c([OH])c([OH])">(mol);
        validate_contains<"N=NC(=S)N">(mol);
        validate_contains<"SC#N">(mol);
        validate_contains<"cC[N+]">(mol);
        validate_contains<"C[O,S;R0][C;R0](=S)">(mol);
        validate_contains<"N[CH2]C#N">(mol);
        validate_contains<"C1(=O)OCC1">(mol);
        validate_contains<"P(=O)([OH])OP(=O)[OH]">(mol);
        validate_contains<"N1CCC1=O">(mol);
        validate_contains<"O=C1[#6]~[#6]C(=O)[#6]~[#6]1">(mol);
        validate_contains<"C=CC=CC=CC=C">(mol);
        validate_contains<"O1CCCCC1OC2CCC3CCCCC3C2">(mol);
        validate_contains<"O=C1NCC2CCCCC21">(mol);
        validate_contains<"O=C1CCCC(N1)=O">(mol);
        validate_contains<"O1CCCCC1C2CCCO2">(mol);
        validate_contains<"[OH]c1cc([OH])cc2OC(C([OH])Cc21)c3cc([OH])c([OH])cc3">(mol);
        validate_contains<"C12OCCC(O1)CC2">(mol);
        validate_contains<"c-N(=O)~O">(mol);
        validate_contains<"[O,N,S][CH2]N1C(=O)CCC1(=O)">(mol);
        validate_contains<"[S,P](=O)(=O)OC">(mol);
        validate_contains<"[CH2]=[CH]-[N,O,S]">(mol);
        validate_contains<"S-C#N">(mol);
        validate_contains<"S=C-[#6,N,O]">(mol);
        validate_contains<"[Cl,Br,I]-N">(mol);
        validate_contains<"C#C-[F,Br,I,Cl]">(mol);
        validate_contains<"C(-[O;H1])(-C#N)">(mol);
        validate_contains<"C(-O)-C-N(=O)=O">(mol);
        validate_contains<"C(=N)-S">(mol);
        validate_contains<"C(=O)-C(-O)-N-C=O">(mol);
        validate_contains<"C(=O)-C(=C)-C(=O)">(mol);
        validate_contains<"C(=O)-C(=N)">(mol);
        validate_contains<"C(=O)-C([F,Br,I,Cl])-C(=O)">(mol);
        validate_contains<"C(=O)-C([F,Br,I,Cl])=C([F,Br,I,Cl])-C(=O)">(mol);
        validate_contains<"C(=O)-C=C-C(=O)">(mol);
        validate_contains<"C(=O)-N(-O)-C(=O)">(mol);
        validate_contains<"C(=O)-N(-S)-C(=O)">(mol);
        validate_contains<"[C;!R](=O)-[N;!R]-[C;!R](=O)">(mol);
        validate_contains<"C(=O)-N-N(=O)">(mol);
        validate_contains<"C(=O)-N-N-C(=O)">(mol);
        validate_contains<"C(=O)-O-N-C(=O)">(mol);
        validate_contains<"C(=O)-S">(mol);
        validate_contains<"C(=O)-S-C(=S)">(mol);
        validate_contains<"C(=S)-S">(mol);
        validate_contains<"[C;H1](=O)">(mol);
        validate_contains<"C([F,Br,I,Cl])=N">(mol);
        validate_contains<"[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]">(mol);
        validate_contains<"C-N=O">(mol);
        validate_contains<"C-S(=O)(=O)-O">(mol);
        validate_contains<"C1-C-C(=O)-O1">(mol);
        validate_contains<"C=C(-S)-S(=O)">(mol);
        validate_contains<"[C;!R]=[C;!R]-[C;!R](-O)">(mol);
        validate_contains<"C=C-C(=N)">(mol);
        validate_contains<"[C;!R]=[C;!R]-[C;!R](=O)">(mol);
        validate_contains<"C=C-C(=O)-C=C">(mol);
        validate_contains<"C=C-C(=S)">(mol);
        validate_contains<"C=C-C(=S)-S">(mol);
        validate_contains<"C=C-C([F,Br,I,Cl])=C([F,Br,I,Cl])">(mol);
        validate_contains<"C=C-C=N">(mol);
        validate_contains<"C=C-N([C,c])([C,c])">(mol);
        validate_contains<"C=C-N(=O)(=O)">(mol);
        validate_contains<"C=C-N(O)(=O)">(mol);
        validate_contains<"C=C-O">(mol);
        validate_contains<"C=C-O-C(=O)">(mol);
        validate_contains<"[C;!R]=[C;!R]-[S;!R]">(mol);
        validate_contains<"C=C-S(=O)(=O)">(mol);
        validate_contains<"C=N-C(-O)">(mol);
        validate_contains<"C=N-C(=O)">(mol);
        validate_contains<"C=N-N">(mol);
        validate_contains<"C=N-N=C">(mol);
        validate_contains<"C=N-O">(mol);
        validate_contains<"C=N-O-C(=O)">(mol);
        validate_contains<"C=N=O">(mol);
        validate_contains<"N#C-S">(mol);
        validate_contains<"N(=O)-C([F,Br,I,Cl])">(mol);
        validate_contains<"[N!H0]-[C!H0]-[N!H0]">(mol);
        validate_contains<"N(O)=C-C-N(=O)">(mol);
        validate_contains<"N-C(-S)=N-S(=O)(=O)">(mol);
        validate_contains<"N-C(=N)-S">(mol);
        validate_contains<"N-C(=O)-O-c1ccccc1">(mol);
        validate_contains<"N-C(=S)">(mol);
        validate_contains<"N-C(=S)-N">(mol);
        validate_contains<"N-C=C">(mol);
        validate_contains<"N-O-C(=O)">(mol);
        validate_contains<"N1-C-C1">(mol);
        validate_contains<"N1-N-C(=O)-N-N1">(mol);
        validate_contains<"N=C(-N)-C(=N)-N">(mol);
        validate_contains<"N=C(-O)-N">(mol);
        validate_contains<"N=C(-S)-C(=N)-N">(mol);
        validate_contains<"N=C([F,Br,I,Cl])">(mol);
        validate_contains<"N=C-C(=O)">(mol);
        validate_contains<"N=C-S">(mol);
        validate_contains<"N=C=N">(mol);
        validate_contains<"N=C=O">(mol);
        validate_contains<"N=C=S">(mol);
        validate_contains<"O-C(=O)-O-N">(mol);
        validate_contains<"O-C([F,Br,I,Cl])=S">(mol);
        validate_contains<"O-C=C">(mol);
        validate_contains<"O-N=C-C=N-O">(mol);
        validate_contains<"P">(mol);
        validate_contains<"S(=O)(=O)-C([F,Br,I,Cl])">(mol);
        validate_contains<"S(=O)(=O)O">(mol);
        validate_contains<"S-C#N">(mol);
        validate_contains<"S-C(=C)-S">(mol);
        validate_contains<"S-C(=N)-N">(mol);
        validate_contains<"S-C(=N)-S">(mol);
        validate_contains<"S-C(=S)-N">(mol);
        validate_contains<"S-N-C(=O)">(mol);
        validate_contains<"[S;H]">(mol);
        validate_contains<"[Si]">(mol);
        validate_contains<"[N;H2]-S">(mol);
        validate_contains<"[N;H]-[C;H]-[N;H]">(mol);
        validate_contains<"c-C(-O)(-O)-c">(mol);
        validate_contains<"c-C(-[O;H])[!O]">(mol);
        validate_contains<"c-S(=O)(=O)-C=C">(mol);
        validate_contains<"c-S(=O)(=O)-O">(mol);
        validate_contains<"n1(=O)c([F,Br,I,Cl])cccc1">(mol);
        validate_contains<"n1c([F,Br,I,Cl])cccc1">(mol);
        validate_contains<"n1c([F,Br,I,Cl])ncccc1">(mol);
        validate_contains<"C1-C-O-C-O-C1">(mol);
        validate_contains<"C=S">(mol);
        validate_contains<"c=S">(mol);
        validate_contains<"S=C-S">(mol);
        validate_contains<"C=C-N(=O)(=O)">(mol);
        validate_contains<"C=C-N(=O)(-O)">(mol);
        validate_contains<"C(=O)-C(=N)">(mol);
        validate_contains<"C=C-S">(mol);
        validate_contains<"C(-S)(-S)(-S)(-S)">(mol);
        validate_contains<"S(-O)(-O)(-O)(-O)">(mol);
        validate_contains<"S(=O)(=O)-C-N(=O)(=O)">(mol);
        validate_contains<"S(=O)(=O)-C-N(=O)(-O)">(mol);
        validate_contains<"N#C-S">(mol);
        validate_contains<"C=N-O">(mol);
        validate_contains<"C=N-S=N">(mol);
        validate_contains<"S-C(=O)-S">(mol);
        validate_contains<"S-C(=N)-S">(mol);
        validate_contains<"S-C(=N)(=N)">(mol);
        validate_contains<"N(=O)(=O)">(mol);
        validate_contains<"N(=O)(-O)">(mol);
        validate_contains<"N(-O)(-O)">(mol);
        validate_contains<"C1-S-C-S1">(mol);
        validate_contains<"N=C=S">(mol);
        validate_contains<"N=C=N">(mol);
        validate_contains<"C1-N=N1">(mol);
        validate_contains<"C([F,Br,Cl,I])-S">(mol);
        validate_contains<"C(-S)(-S)=C(-S)(-S)">(mol);
        validate_contains<"C(=N)-C=N-O">(mol);
        validate_contains<"C(-C#N)(-C#N)">(mol);
        validate_contains<"C-O-N=O">(mol);
        validate_contains<"C=C-N">(mol);
        validate_contains<"C=N-S(=O)(=O)">(mol);
        validate_contains<"N-N=O">(mol);
        validate_contains<"N-C(=O)-S">(mol);
        validate_contains<"S-C#N">(mol);
        validate_contains<"N-C(=O)-C(=O)-N">(mol);
        validate_contains<"c=N">(mol);
        validate_contains<"C1-O-C-O-C1">(mol);
        validate_contains<"C=N-N-S(=O)(=O)">(mol);
        validate_contains<"[C;!R](=O)-[C;!R](=O)">(mol);
        validate_contains<"C(=O)-N-O-C(=O)">(mol);
        validate_contains<"C1(=O)-C=C-C(=O)-C=C1">(mol);
        validate_contains<"B">(mol);
        validate_contains<"ONO">(mol);
        validate_contains<"ON(~O)~O">(mol);
        validate_contains<"C1(O)CC(O)CO1">(mol);
        validate_contains<"C1(O)C(O)CCO1">(mol);
        validate_contains<"C1(O)CC(O)CCO1">(mol);
        validate_contains<"N=[N+]=N">(mol);
        validate_contains<"N-C#N">(mol);
        validate_contains<"[C;!R]=N-N=[C;!R]">(mol);
        validate_contains<"N(C)(C)-[C;H2]-[C;H2]([F,Br,I,Cl])">(mol);
        validate_contains<"[C;a]-[N;H1]-[N;H2]">(mol);
        validate_contains<"[C;a]-[C;H2]([F,Br,I,Cl])">(mol);
        validate_contains<"[S;!R]-[C;!R]-[O;!R]">(mol);
        validate_contains<"N1-C-O1">(mol);
        validate_contains<"[C;H2]-O-S(=O)(=O)-C">(mol);
        validate_contains<"C1-C-O1">(mol);
        validate_contains<"c1-[N;a]-c([F,Br,I,Cl])-c-c-c1">(mol);
        validate_contains<"C=[C;H0]([F,Br,I,Cl])([F,Br,I,Cl])">(mol);
        validate_contains<"[N+]#N">(mol);
        validate_contains<"ON#C">(mol);
        validate_contains<"OOO">(mol);
        validate_contains<"OO">(mol);
        validate_contains<"Cl(~O)(~O)(~O)">(mol);
        validate_contains<"Cl(~O)(~O)(~O)(~O)">(mol);
        validate_contains<"C(~O)(~O)(~[OH])">(mol);
        validate_contains<"N[F,I,Br,Cl]">(mol);
        validate_contains<"P[F,I,Br,Cl]">(mol);
        validate_contains<"S[F,I,Br,Cl]">(mol);
        validate_contains<"C=C-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"C=C-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"C=C-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"C=C-[C;H2]([F,Br,I,Cl])">(mol);
        validate_contains<"CN=C=O">(mol);
        validate_contains<"CN=C=S">(mol);
        validate_contains<"C(=O)([F,I,Br,Cl])">(mol);
        validate_contains<"[N;!R]=[N;!R]">(mol);
        validate_contains<"[S;!R]-[S;!R]">(mol);
        validate_contains<"C1NC(=O)OC(=O)C1">(mol);
        validate_contains<"[C;!R]C(=N)S[C;!R]">(mol);
        validate_contains<"[C;!R]C(=N)O[C;!R]">(mol);
        validate_contains<"C=C=C">(mol);
        validate_contains<"C(=O)-O-C(=O)-[!N]">(mol);
        validate_contains<"C=C-C#N">(mol);
        validate_contains<"[C;H2]=C-C(=O)">(mol);
        validate_contains<"[C;D1&H3,D2&H2,D3&H1,D4]-[C;H1]=C-C(=O)-C">(mol);
        validate_contains<"[C;D1&H3,D2&H2,D3&H1,D4]-[C;H1]=C-C(=O)-O-C-C">(mol);
        validate_contains<"C#C-C#N">(mol);
        validate_contains<"C#C-C(=O)">(mol);
        validate_contains<"C#C-C(=O)-O-C-C">(mol);
        validate_contains<"C#C-[C;!H0](=O)">(mol);
        validate_contains<"C#C-C(=O)[N,!H0]">(mol);
        validate_contains<"C=C-[C;H1](=O)">(mol);
        validate_contains<"[O;D2]-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"[C;!H3]-O-S(=O)(=O)O">(mol);
        validate_contains<"[C;!a]-S(=O)(=O)[O-]">(mol);
        validate_contains<"[N;!a]-S(=O)(=O)[O-]">(mol);
        validate_contains<"C(=O)-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"C(=O)-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"C(=O)-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"C(=O)-[C;H2]([F,Br,I,Cl])">(mol);
        validate_contains<"S(=O)(=O)-[C;!H0]([F,Br,I,Cl])([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"S(=O)(=O)-[C;!H0]([F,Br,I,Cl])([!F,!Br,!I,!Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"S(=O)(=O)-[C;H1]([F,Br,I,Cl])([!F,!Br,!I,!Cl])">(mol);
        validate_contains<"S(=O)(=O)-[C;H2]([F,Br,I,Cl])">(mol);
        validate_contains<"[O,N,S]1CC1">(mol);
        validate_contains<"O=COC=O">(mol);
        validate_contains<"N=C=O">(mol);
        validate_contains<"S=C=N">(mol);
        validate_contains<"O~O">(mol);
        validate_contains<"[Si]~N">(mol);
        validate_contains<"P[S,N]">(mol);
        validate_contains<"[N;R0]=[N;R0]=[C;R0]">(mol);
        validate_contains<"[N+]#N-*">(mol);
        validate_contains<"[C,c]-[S;D2]-[O,N]">(mol);
        validate_contains<"[Cl,Br,I,F][S,P,Si,N]">(mol);
        validate_contains<"[Cl,Br,I]CC[S,N]">(mol);
        validate_contains<"[Si]-O-*">(mol);
        validate_contains<"S(O[C,c,N,n])(~O)[C,c,N,n]">(mol);
        validate_contains<"[N;R0](~N)~O">(mol);
        validate_contains<"N(~O)(~O)(~O)-*">(mol);
        validate_contains<"[N+]([O-])(=C)-*">(mol);
        validate_contains<"N([H,C;X4])([H,C;X4])-[C;R0;X4]-N[C;X4]([H,C;X4])([H,C;X4])">(mol);
        validate_contains<"[C;R0](=[C;R0])-[OH]">(mol);
        validate_contains<"[Be,B,Al,Ti,Cr,Mn,Fe,Co,Ni,Cu,Pd,Ag,Sn,Pt,Au,Hg,Pb,Bi,As,Sb,Gd,Te]">(mol);
        validate_contains<"O=CON1C(=O)CCC1=O">(mol);
        validate_contains<"O=COn1cncc1">(mol);
        validate_contains<"Fc1c(OC=O)c(F)c(F)c(F)c1F">(mol);
        validate_contains<"[S;D2][C;R0](C)(C)C">(mol);
        validate_contains<"[CH]=O">(mol);
        validate_contains<"[S;D2;R0]-[S;D2]">(mol);
        validate_contains<"[SH]">(mol);
        validate_contains<"[S,s;D2]C[S,s;D2]">(mol);
        validate_contains<"[n+]-C">(mol);
        validate_contains<"[NH]=C([NH2])c">(mol);
        validate_contains<"O=CN[OH]">(mol);
        validate_contains<"[NH;R0][NH;R0]">(mol);
        validate_contains<"P(=O)(O[H,C])O[H,C]">(mol);
        validate_contains<"O=CC([Cl,Br,I,F])([Cl,Br,I,F])[Cl,Br,I,F]">(mol);
        validate_contains<"[#6]SC(=O)">(mol);
        validate_contains<"C(=O)Cl">(mol);
        validate_contains<"C(=O)[Cl,Br]">(mol);
        validate_contains<"[CX4][OH]">(mol);
        validate_contains<"[CX4;H2][OH]">(mol);
        validate_contains<"[CX4;H][OH]">(mol);
        validate_contains<"[CX4;H0][OH]">(mol);
        validate_contains<"[#6][CH]=O">(mol);
        validate_contains<"[CX4][Cl,Br,I]">(mol);
        validate_contains<"C#C">(mol);
        validate_contains<"[Cl,Br,I][CH]C=C">(mol);
        validate_contains<"cCC(=O)[OH]">(mol);
        validate_contains<"C([Cl,Br,I])([Cl,Br,I])([Cl,Br,I])C(=O)[OH]">(mol);
        validate_contains<"[#6]C(=O)C(=O)[OH]">(mol);
        validate_contains<"C=CC(=O)[OH]">(mol);
        validate_contains<"N#C-C(=O)[OH]">(mol);
        validate_contains<"N(~O)(~O)-C-C(=O)[OH]">(mol);
        validate_contains<"[OH][CX4][Cl,Br,I]">(mol);
        validate_contains<"O(C(=O)[#6])(C(=O)[#6])">(mol);
        validate_contains<"c[NH,NH2]">(mol);
        validate_contains<"c-[ND3]([#6])[#6]">(mol);
        validate_contains<"c">(mol);
        validate_contains<"c-[F,Cl,Br,I]">(mol);
        validate_contains<"c-[Br,I]">(mol);
        validate_contains<"cB([OH])[OH]">(mol);
        validate_contains<"[N;R0;D2]~[N;R0]~[N;R0;D1]">(mol);
        validate_contains<"[OH][CX4][CX4][Cl,Br,I]">(mol);
        validate_contains<"[CX4]([OH])[CX4]C(=O)[OH]">(mol);
        validate_contains<"[CX4]([Cl,Br,I])[CX4]C(=O)[OH]">(mol);
        validate_contains<"[#6]C(=O)[CX4]C(=O)[OH]">(mol);        
        validate_contains<"C(=O)[OH]">(mol);
        validate_contains<"O([CX4,c])C(=O)O[CX4,c]">(mol);
        validate_contains<"ClC(=O)O[CX4,c]">(mol);
        validate_contains<"[CX4,c]-C#N">(mol);
        validate_contains<"N[CX4]C(=O)N[CX4]C(=O)">(mol);
        validate_contains<"[S;D2]-[S;D2]">(mol);
        validate_contains<"O1[CX4][CX4]1">(mol);
        validate_contains<"[OH][CH,CH2]O[CX4,c]">(mol);
        validate_contains<"O([#6])-C([#6])([#6])-[OH]">(mol);
        validate_contains<"C=N[NH2]">(mol);
        validate_contains<"C=NC=O">(mol);
        validate_contains<"C#N-[#6]">(mol);
        validate_contains<"O=C=N-[#6]">(mol);
        validate_contains<"S=C=N-[#6]">(mol);
        validate_contains<"O([CX4,c])-C([CX4,c])([CX4,c])-O([CX4,c])">(mol);
        validate_contains<"[CX4,c]C(=O)[CX4,c]">(mol);
        validate_contains<"OC(=O)CC(=O)[OH]">(mol);
        validate_contains<"[r8,r9,r10,r11,r12,r13,r14]">(mol);
        validate_contains<"*-C#N">(mol);
        validate_contains<"*-N(=O)(~O)">(mol);
        validate_contains<"[N;D4]">(mol);
        validate_contains<"OC(=O)C(=O)O">(mol);
        validate_contains<"C=[N;R0]-[OH]">(mol);
        validate_contains<"O~O">(mol);
        validate_contains<"c[OH]">(mol);
        validate_contains<"[PX3](=O)(~O)~[OH]">(mol);
        validate_contains<"[PX3](=O)(~O)~O-[#6]">(mol);
        validate_contains<"[PX4](=O)(~O)(~O)~[OH]">(mol);
        validate_contains<"[PX4](=O)(~O)(~O)~O-[#6]">(mol);
        validate_contains<"[PX4](=O)(~O)~O">(mol);
        validate_contains<"O=P(~O)(~O)(~O)">(mol);
        validate_contains<"[CX4][NH2]">(mol);
        validate_contains<"C14~*~*~*~*~C~1~*~*~C2~C3~*~*~*~C~3~*~*~C~2~4">(mol);
        validate_contains<"[#6][SD3](~O)[#6]">(mol);
        validate_contains<"[#6][SD4](~O)(~O)[#6]">(mol);
        validate_contains<"[#6][SD4](~O)(~O)N">(mol);
        validate_contains<"[#6]S(~O)(~O)[OH]">(mol);
        validate_contains<"[#6]S(~O)(~O)O[#6]">(mol);
        validate_contains<"[#6][SD4](~O)(~O)[Cl,Br]">(mol);
        validate_contains<"[ND3]([CX4])([CX4])[CX4]">(mol);
        validate_contains<"s1c(=N)nnc1[S,N]">(mol);
        validate_contains<"[CX4][SH]">(mol);
        validate_contains<"c[SH]">(mol);
        validate_contains<"NC(=S)N">(mol);
        validate_contains<"[D2R0]-[D2R0]-[D2R0]-[D2R0]-[D2R0]">(mol);
        validate_contains<"NC(=O)N">(mol);


        // ! in bonds -> bug
        //validate_contains<"O=C-N!@N">(mol);
        //validate_contains<"C!@N=*">(mol);
        //validate_contains<"[C;R0;X4]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH]!@[CH][CH2]">(mol);
        //validate_contains<"[C;R0;X4]!@[CH2]!@[CH2]!@[CH2]!@[CH2]!@[CH2]!@[C;R0;X4]">(mol);
        //validate_contains<"[C;R0;X4]!@[CX4]!@[CX4]!@[CX4]!@[CX4]!@[CX4]!@[C;R0;X4]">(mol);
        //validate_contains<"[c,C]1(~[O;D1])~*!-*~[c,C](~[O;D1])~*!-*~1">(mol);

        // recursive SMARTS -> not yet supported
        //validate_contains<"[#6]C(=O)[CX4]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol);
        //validate_contains<"O=CC([$(C(F)(F)F),$(C#N),$(Cl)])=C">(mol);
        //validate_contains<"[$(N!@N),$([N;R0]=[N;R0])]">(mol);
        //validate_contains<"[$(N!@S),$(N~C(~S)-N)]">(mol);
        //validate_contains<"[$(C-[Cl,Br,I]),$(O=C-[Cl,Br,I,F]),$(O=C([CH,CH2][Cl,Br,I,F])[O,C]),$(C~O~[Cl,Br,I,F][CH,CH2]),$(n1c([Cl,Br,I,F])nccc1);!$(C=C-[Cl,Br,I]);!$(ClC-[Cl,Br,I,F])]">(mol);
        //validate_contains<"[!$([C,c]-N(=O)~O);$([!O]~[N;R0]=O)]">(mol);
        //validate_contains<"N#C[C;R0;X4]O[!$(O=[C,S])]">(mol);
        //validate_contains<"[C;R0](=[C;R0])-[S,O,N;R0][!$(O=[C,S])]">(mol);
        //validate_contains<"[$([S,s]~[S,s]~[C,c]=S),$([S,s]~[C,c](=S)~[S,s,N]),$([S;D2;R0]-S~O)]">(mol);
        //validate_contains<"[$(O=C[CH](C=O)C=O),$(N#C[CH](-C=O)-C=O)]">(mol);
        //validate_contains<"[$(N#C-C=[CH][C,c]),$([CH](=[C;R0]-[CH]=O)),$([CH](=[C,R]-C(=O)-C));!$([CH]1=CC(=O)C=CC1=*);!$([CH]1=CC(=O)C(=[N,O])C=C1);!$([CH](=C-C=O)-C=O)]">(mol);
        //validate_contains<"[$(N#C-C#C[C,c]),$(C#C-[CH]=O),$(C(#C-C(=O)-[C,c]))]">(mol);
        //validate_contains<"[$(N#CSc1sc(nc1)N),$([S,Se]1C(N)C(=O)[#6][#6]1)]">(mol);
        //validate_contains<"[!$(O=[C,S])][N;R0]=[C;R0]([C,c])[C,c]">(mol);
        //validate_contains<"[$(O([CX4,c])!@[CH,CH2]!@O[CX4,c])]">(mol);
        //validate_contains<"[$(C#N),$([C,N,S]=O)][CH2,CH][$([C,N,S]=O),$(C#N)]">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]">(mol);validate_contains<"[O-]-[C;R0]=[C;R0]">(mol);
        //validate_contains<"[ND3]([CX4,c,H])([CX4,c,H])[CX4][$([CH]),$(C([CX4,c]))]=O">(mol);
        //validate_contains<"[OH][CX4][$([CH]),$(C([CX4,c]))]=O">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C-[$([CH]),$(C([CX4,c]))]=O">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C-C#N">(mol);
        //validate_contains<"[Cl,Br,I][$([CX4][CH]=O),$([CX4]C(=O)[CX4,c])]">(mol);
        //validate_contains<"[OH][CX4][$([NH2]),$([NH][CX4]),$(N([CX4])[CX4])]">(mol);
        //validate_contains<"O=C([c,CX4])[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol);
        //validate_contains<"[$([NH2][CX4]),$([$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]),$([ND3]([CX4])([CX4])[CX4])]">(mol);
        //validate_contains<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4]);!$(NC=O)][CX4]C(=O)[OH]">(mol);
        //validate_contains<"Fcaa[F,Cl,Br,I,$([C,N,S]=O)]">(mol);
        //validate_contains<"[N;D2]([C,c;!$(C=[O,S,N])])=[N;D2]-[C,c;!$(C=[O,S,N])]">(mol);
        //validate_contains<"[OH][CX4][CX4][$([NH2]),$([NH][CX4]),$(N([CX4])[CX4])]">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=C[CX4]C(=O)[OH]">(mol);
        //validate_contains<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol);
        //validate_contains<"[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]C(=O)[OH]">(mol);
        //validate_contains<"[$([CX4,c][CH]=O),$([CX4,c]C(=O)[CX4,c])]">(mol);
        //validate_contains<"[$([S;D2]([CX4,c])!@[CH,CH2]!@[S;D2][CX4,c])]">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3][$([NH2]),$([NH][CX4]),$([N;R0]([CX4])[CX4])]);!$(CC=[O,S,N]);!$(C[O,S])]">(mol);
        //validate_contains<"[$([CX3]C(=O)[CX4,c]);!$(CC=[S,N]);!$(C[O,S,N])]=[$([CX3]C(=O)[CX4,c]);!$(CC=[S,N]);!$(C[N,S])]">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3][OH]);!$(CC=[O,S,N]);!$(C[N,S])]">(mol);
        //validate_contains<"[$([CX3]);!$(CC=[O,S,N]);!$(C[O,S,N])]=[$([CX3]O[CX4]);!$(CC=[O,S,N]);!$(C[N,S])]">(mol);
        //validate_contains<"[$([#6]);!$(C=[O,S,N])]C(=O)O[$([#6]);!$(C=[O,S,N])]">(mol);
        //validate_contains<"[$(O([$([CX4,c]);!$(C[O,N,S])])[$([CX4,c]);!$(C[O,N,S])]);!$(O1[CX4][CX4]1)]">(mol);
        //validate_contains<"[OH][$([NX3]([C;!$(C=[O,S,N])])[C;!$(C=[O,S,N])]),$([NH][CX4])]">(mol);
        //validate_contains<"[$([NH;R0]([C;!$(C=[O,S,N])]))][$([NH;R0][C;!$(C=[O,S,N])])]">(mol);
        //validate_contains<"[$([C;R1]);!$(C(N)N)](=O)@[$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol);
        //validate_contains<"[$([C;R1]);!$(C(O)N);!$(C(O)O)](=O)@[$(O);!$(O(C=O))]">(mol);
        //validate_contains<"[$([C;R0]=[N;R0]);!$(C(~[N,O,S])(~[N,O,S]));!$([C;R0]=[N;R0]~[N,O,n])]">(mol);
        //validate_contains<"[$([NH]([CX4])[CX4]);!$([NH]([CX4])[CX4][O,N]);!$([NH]([CX4])[CX4][O,N])]">(mol);
        //validate_contains<"S=C([c,CX4])[$([NH2]),$([NH][c,CX4]),$(N([c,CX4])[c,CX4])]">(mol);
        //validate_contains<"[$([#6]);!$(C=[O,S,N])]C(=S)O[$([#6]);!$(C=[O,S,N])]">(mol);
        //validate_contains<"[$([#16;D2]([$([CX4,c]);!$(C[O,N,S])])[$([CX4,c]);!$(C[O,N,S])]);!$([#16;D2]1[CX4][CX4]1);!$(s1aaaa1)]">(mol);
        //validate_contains<"[$([CX4,c][CH]=S),$([CX4,c]C(=S)[CX4,c])]">(mol);
    }
}
