#include <Kitimar/CTSmarts/CTSmarts.hpp>



#ifdef KITIMAR_WITH_OPENBABEL
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <openbabel/parsmart.h>
#endif

#include <gtest/gtest.h>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

template <ctll::fixed_string SMARTS, typename AtomExpr/*, int AtomIndex = 0*/>
constexpr void test_atom_expr()
{
    constexpr auto AtomIndex = 0;
    auto smarts = Smarts<SMARTS>{};
    static_assert(ctll::size(smarts.atoms) > AtomIndex);
    auto expr = get<AtomIndex>(smarts.atoms);
    //static_assert(std::is_same_v<decltype(expr), AtomExpr>);
    //static_assert(std::is_same_v<decltype(expr), decltype(optimizeExpression(AtomExpr{}))>);
    static_assert(std::is_same_v<decltype(expr), Atom<AtomIndex, AtomExpr>>);


}

template <ctll::fixed_string SMARTS, typename BondExpr, int BondIndex = 0>
constexpr void test_bond_expr()
{
    auto smarts = Smarts<SMARTS>{};
    static_assert(ctll::size(smarts.bonds) > BondIndex);
    auto expr = get<BondIndex>(smarts.bonds);
    static_assert(std::is_same_v<decltype(expr), Bond<BondIndex, 0, 1, BondExpr>>);
}


template <ctll::fixed_string SMARTS, typename BondT, int BondIndex = 0>
constexpr void test_bond()
{
    auto smarts = Smarts<SMARTS>{};
    static_assert(ctll::size(smarts.bonds) > BondIndex);
    auto expr = get<BondIndex>(smarts.bonds);
    static_assert(std::is_same_v<decltype(expr), BondT>);
}


template <ctll::fixed_string SMARTS, typename Error>
constexpr void test_parse_error()
{
    auto smarts = Smarts<SMARTS, true>();
    static_assert(!smarts.valid);
    auto error = smarts.error;
    static_assert(std::is_same_v<decltype(error), Error>);
}

template <ctll::fixed_string SMARTS, typename Error>
constexpr void test_error()
{
    auto smarts = Smarts<SMARTS, true>();
    static_assert(smarts.valid);
    auto error = smarts.error;
    static_assert(std::is_same_v<decltype(error), Error>);
}







TEST(TestCTSmarts, AliphaticAtom)
{
    // 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
    test_atom_expr<"B",  AliphaticAtom<5>>();
    test_atom_expr<"C",  AliphaticAtom<6>>();
    test_atom_expr<"N",  AliphaticAtom<7>>();
    test_atom_expr<"O",  AliphaticAtom<8>>();
    test_atom_expr<"P",  AliphaticAtom<15>>();
    test_atom_expr<"S",  AliphaticAtom<16>>();
    test_atom_expr<"F",  AliphaticAtom<9>>();
    test_atom_expr<"Cl", AliphaticAtom<17>>();
    test_atom_expr<"Br", AliphaticAtom<35>>();
    test_atom_expr<"I",  AliphaticAtom<53>>();
}

TEST(TestCTSmarts, AromaticAtom)
{
    // 'b' | 'c' | 'n' | 'o' | 's' | 'p'
    test_atom_expr<"b",  AromaticAtom<5>>();
    test_atom_expr<"c",  AromaticAtom<6>>();
    test_atom_expr<"n",  AromaticAtom<7>>();
    test_atom_expr<"o",  AromaticAtom<8>>();
    test_atom_expr<"p",  AromaticAtom<15>>();
    test_atom_expr<"s",  AromaticAtom<16>>();
}

TEST(TestCTSmarts, AnyAtom)
{
    test_atom_expr<"*",  AnyAtom>();
    test_atom_expr<"a",  AnyAromatic>();
    test_atom_expr<"A",  AnyAliphatic>();
}

TEST(TestCTSmarts, BracketAliphaticAtom)
{
    test_atom_expr<"[H]",  AliphaticAtom<1> >();
    test_atom_expr<"[He]", AliphaticAtom<2> >();
    test_atom_expr<"[Li]", AliphaticAtom<3> >();
    test_atom_expr<"[Be]", AliphaticAtom<4> >();
    test_atom_expr<"[B]",  AliphaticAtom<5> >();
    test_atom_expr<"[C]",  AliphaticAtom<6> >();
    test_atom_expr<"[N]",  AliphaticAtom<7> >();
    test_atom_expr<"[O]",  AliphaticAtom<8> >();
    test_atom_expr<"[F]",  AliphaticAtom<9> >();
    test_atom_expr<"[Ne]", AliphaticAtom<10> >();
    test_atom_expr<"[Na]", AliphaticAtom<11> >();
    test_atom_expr<"[Mg]", AliphaticAtom<12> >();
    test_atom_expr<"[Al]", AliphaticAtom<13> >();
    test_atom_expr<"[Si]", AliphaticAtom<14> >();
    test_atom_expr<"[P]",  AliphaticAtom<15> >();
    test_atom_expr<"[S]",  AliphaticAtom<16> >();
    test_atom_expr<"[Cl]", AliphaticAtom<17> >();
    test_atom_expr<"[Ar]", AliphaticAtom<18> >();
    test_atom_expr<"[K]",  AliphaticAtom<19> >();
    test_atom_expr<"[Ca]", AliphaticAtom<20> >();
    test_atom_expr<"[Sc]", AliphaticAtom<21> >();
    test_atom_expr<"[Ti]", AliphaticAtom<22> >();
    test_atom_expr<"[V]",  AliphaticAtom<23> >();
    test_atom_expr<"[Cr]", AliphaticAtom<24> >();
    test_atom_expr<"[Mn]", AliphaticAtom<25> >();
    test_atom_expr<"[Fe]", AliphaticAtom<26> >();
    test_atom_expr<"[Co]", AliphaticAtom<27> >();
    test_atom_expr<"[Ni]", AliphaticAtom<28> >();
    test_atom_expr<"[Cu]", AliphaticAtom<29> >();
    test_atom_expr<"[Zn]", AliphaticAtom<30> >();
    test_atom_expr<"[Ga]", AliphaticAtom<31> >();
    test_atom_expr<"[Ge]", AliphaticAtom<32> >();
    test_atom_expr<"[As]", AliphaticAtom<33> >();
    test_atom_expr<"[Se]", AliphaticAtom<34> >();
    test_atom_expr<"[Br]", AliphaticAtom<35> >();
    test_atom_expr<"[Kr]", AliphaticAtom<36> >();
    test_atom_expr<"[Rb]", AliphaticAtom<37> >();
    test_atom_expr<"[Sr]", AliphaticAtom<38> >();
    test_atom_expr<"[Y]",  AliphaticAtom<39> >();
    test_atom_expr<"[Zr]", AliphaticAtom<40> >();
    test_atom_expr<"[Nb]", AliphaticAtom<41> >();
    test_atom_expr<"[Mo]", AliphaticAtom<42> >();
    test_atom_expr<"[Tc]", AliphaticAtom<43> >();
    test_atom_expr<"[Ru]", AliphaticAtom<44> >();
    test_atom_expr<"[Rh]", AliphaticAtom<45> >();
    test_atom_expr<"[Pd]", AliphaticAtom<46> >();
    test_atom_expr<"[Ag]", AliphaticAtom<47> >();
    test_atom_expr<"[Cd]", AliphaticAtom<48> >();
    test_atom_expr<"[In]", AliphaticAtom<49> >();
    test_atom_expr<"[Sn]", AliphaticAtom<50> >();
    test_atom_expr<"[Sb]", AliphaticAtom<51> >();
    test_atom_expr<"[Te]", AliphaticAtom<52> >();
    test_atom_expr<"[I]",  AliphaticAtom<53> >();
    test_atom_expr<"[Xe]", AliphaticAtom<54> >();
    test_atom_expr<"[Cs]", AliphaticAtom<55> >();
    test_atom_expr<"[Ba]", AliphaticAtom<56> >();
    test_atom_expr<"[La]", AliphaticAtom<57> >();
    test_atom_expr<"[Ce]", AliphaticAtom<58> >();
    test_atom_expr<"[Pr]", AliphaticAtom<59> >();
    test_atom_expr<"[Nd]", AliphaticAtom<60> >();
    test_atom_expr<"[Pm]", AliphaticAtom<61> >();
    test_atom_expr<"[Sm]", AliphaticAtom<62> >();
    test_atom_expr<"[Eu]", AliphaticAtom<63> >();
    test_atom_expr<"[Gd]", AliphaticAtom<64> >();
    test_atom_expr<"[Tb]", AliphaticAtom<65> >();
    test_atom_expr<"[Dy]", AliphaticAtom<66> >();
    test_atom_expr<"[Ho]", AliphaticAtom<67> >();
    test_atom_expr<"[Er]", AliphaticAtom<68> >();
    test_atom_expr<"[Tm]", AliphaticAtom<69> >();
    test_atom_expr<"[Yb]", AliphaticAtom<70> >();
    test_atom_expr<"[Lu]", AliphaticAtom<71> >();
    test_atom_expr<"[Hf]", AliphaticAtom<72> >();
    test_atom_expr<"[Ta]", AliphaticAtom<73> >();
    test_atom_expr<"[W]",  AliphaticAtom<74> >();
    test_atom_expr<"[Re]", AliphaticAtom<75> >();
    test_atom_expr<"[Os]", AliphaticAtom<76> >();
    test_atom_expr<"[Ir]", AliphaticAtom<77> >();
    test_atom_expr<"[Pt]", AliphaticAtom<78> >();
    test_atom_expr<"[Au]", AliphaticAtom<79> >();
    test_atom_expr<"[Hg]", AliphaticAtom<80> >();
    test_atom_expr<"[Tl]", AliphaticAtom<81> >();
    test_atom_expr<"[Pb]", AliphaticAtom<82> >();
    test_atom_expr<"[Bi]", AliphaticAtom<83> >();
    test_atom_expr<"[Po]", AliphaticAtom<84> >();
    test_atom_expr<"[At]", AliphaticAtom<85> >();
    test_atom_expr<"[Rn]", AliphaticAtom<86> >();
    test_atom_expr<"[Fr]", AliphaticAtom<87> >();
    test_atom_expr<"[Ra]", AliphaticAtom<88> >();
    test_atom_expr<"[Ac]", AliphaticAtom<89> >();
    test_atom_expr<"[Th]", AliphaticAtom<90> >();
    test_atom_expr<"[Pa]", AliphaticAtom<91> >();
    test_atom_expr<"[U]",  AliphaticAtom<92> >();
    test_atom_expr<"[Np]", AliphaticAtom<93> >();
    test_atom_expr<"[Pu]", AliphaticAtom<94> >();
    test_atom_expr<"[Am]", AliphaticAtom<95> >();
    test_atom_expr<"[Cm]", AliphaticAtom<96> >();
    test_atom_expr<"[Bk]", AliphaticAtom<97> >();
    test_atom_expr<"[Cf]", AliphaticAtom<98> >();
    test_atom_expr<"[Es]", AliphaticAtom<99> >();
    test_atom_expr<"[Fm]", AliphaticAtom<100> >();
    test_atom_expr<"[Md]", AliphaticAtom<101> >();
    test_atom_expr<"[No]", AliphaticAtom<102> >();
    test_atom_expr<"[Lr]", AliphaticAtom<103> >();
}

TEST(TestCTSmarts, BracketAromaticAtom)
{
    test_atom_expr<"[b]",  AromaticAtom<5> >();
    test_atom_expr<"[c]",  AromaticAtom<6> >();
    test_atom_expr<"[n]",  AromaticAtom<7> >();
    test_atom_expr<"[o]",  AromaticAtom<8> >();
    test_atom_expr<"[p]",  AromaticAtom<15> >();
    test_atom_expr<"[s]",  AromaticAtom<16> >();
    test_atom_expr<"[as]", AromaticAtom<33> >();
    test_atom_expr<"[se]", AromaticAtom<34> >();
}

TEST(TestCTSmarts, BracketAnyAtom)
{
    test_atom_expr<"[*]", AnyAtom >();
    test_atom_expr<"[a]", AnyAromatic >();
    test_atom_expr<"[A]", AnyAliphatic >();
}

TEST(TestCTSmarts, Isotope)
{
    test_atom_expr<"[1]",   Isotope<1> >();
    test_atom_expr<"[2]",   Isotope<2> >();
    test_atom_expr<"[0]",   Isotope<0> >();
    test_atom_expr<"[42]",  Isotope<42> >();
    test_atom_expr<"[123]", Isotope<123> >();
}

TEST(TestCTSmarts, Element)
{
    test_atom_expr<"[#1]",   Element<1> >();
    test_atom_expr<"[#2]",   Element<2> >();
    test_atom_expr<"[#0]",   Element<0> >();
    test_atom_expr<"[#42]",  Element<42> >();
    test_atom_expr<"[#123]", Element<123> >();
}

TEST(TestCTSmarts, Degree)
{
    test_atom_expr<"[D]",   Degree<1> >();
    test_atom_expr<"[D1]",  Degree<1> >();
    test_atom_expr<"[D2]",  Degree<2> >();
    test_atom_expr<"[D0]",  Degree<0> >();
    test_atom_expr<"[D42]", Degree<42> >();
}

TEST(TestCTSmarts, Valence)
{
    test_atom_expr<"[v]",   Valence<1> >();
    test_atom_expr<"[v1]",  Valence<1> >();
    test_atom_expr<"[v2]",  Valence<2> >();
    test_atom_expr<"[v0]",  Valence<0> >();
    test_atom_expr<"[v42]", Valence<42> >();
}

TEST(TestCTSmarts, Connectivity)
{
    test_atom_expr<"[X]",   Connectivity<1> >();
    test_atom_expr<"[X1]",  Connectivity<1> >();
    test_atom_expr<"[X2]",  Connectivity<2> >();
    test_atom_expr<"[X0]",  Connectivity<0> >();
    test_atom_expr<"[X42]", Connectivity<42> >();
}

TEST(TestCTSmarts, TotalH)
{
    test_atom_expr<"[H]",   AliphaticAtom<1> >();
    test_atom_expr<"[rH]",  And<Cyclic, AliphaticAtom<1>> >();
    test_atom_expr<"[H0]", TotalH<0> >();
    test_atom_expr<"[H1]", TotalH<1> >();
    test_atom_expr<"[H3]", TotalH<3> >();
    test_atom_expr<"[*H]",  And<AnyAtom, TotalH<1>> >();
    test_atom_expr<"[*H0]", And<AnyAtom, TotalH<0>> >();
    test_atom_expr<"[*H1]", And<AnyAtom, TotalH<1>> >();
    test_atom_expr<"[*H3]", And<AnyAtom, TotalH<3>> >();
    test_atom_expr<"[*,H]",  Or<AnyAtom, TotalH<1>> >();
    test_atom_expr<"[*,H0]", Or<AnyAtom, TotalH<0>> >();
    test_atom_expr<"[*,H1]", Or<AnyAtom, TotalH<1>> >();
    test_atom_expr<"[*,H3]", Or<AnyAtom, TotalH<3>> >();

    static_assert(!SmartsActions::isTotalHExpr(AliphaticAtom<1>{}));
    static_assert(!SmartsActions::isTotalHExpr(Cyclic{}));
    static_assert(SmartsActions::isTotalHExpr(AnyAtom{}));
    static_assert(SmartsActions::isTotalHExpr(AnyAliphatic{}));
    static_assert(SmartsActions::isTotalHExpr(AnyAromatic{}));
    static_assert(SmartsActions::isTotalHExpr(AliphaticAtom<6>{}));
    static_assert(SmartsActions::isTotalHExpr(AromaticAtom<6>{}));
    static_assert(SmartsActions::isTotalHExpr(Or<AliphaticAtom<1>, AliphaticAtom<6>>{}));
    static_assert(SmartsActions::isTotalHExpr(And<AnyAtom, AliphaticAtom<1>>{}));
    static_assert(SmartsActions::isTotalHExpr(Or<And<AliphaticAtom<7>, AliphaticAtom<1>>, And<AliphaticAtom<7>, TotalH<2>>>{}));

    test_atom_expr<"[H,R]", Or<AliphaticAtom<1>, Cyclic> >();
    test_atom_expr<"[C,H]", Or<AliphaticAtom<6>, TotalH<1>> >();
    test_atom_expr<"[H,C]", Or<TotalH<1>, AliphaticAtom<6>> >();
    test_atom_expr<"[NH,NH2]", Or<And<AliphaticAtom<7>, TotalH<1>>, And<AliphaticAtom<7>, TotalH<2>>> >();
}

TEST(TestCTSmarts, ImplicitH)
{
    test_atom_expr<"[h]",  ImplicitH<1> >();
    test_atom_expr<"[h0]", ImplicitH<0> >();
    test_atom_expr<"[h1]", ImplicitH<1> >();
    test_atom_expr<"[h3]", ImplicitH<3> >();
}

TEST(TestCTSmarts, CyclicAtom)
{
    test_atom_expr<"[R]", Cyclic >();
    test_atom_expr<"[r]", Cyclic >();
    test_atom_expr<"[x]", Cyclic >();
}

TEST(TestCTSmarts, AcyclicAtom)
{
    //test_atom_expr<"[R0]", Acyclic >();  FIXME
    test_atom_expr<"[r0]", Acyclic >();
    test_atom_expr<"[x0]", Acyclic >();
}

TEST(TestCTSmarts, RingCount)
{
    test_atom_expr<"[R1]",  RingCount<1> >();
    test_atom_expr<"[R2]",  RingCount<2> >();
    test_atom_expr<"[R42]", RingCount<42> >();
}

TEST(TestCTSmarts, RingSize)
{
    test_atom_expr<"[r1]",  RingSize<1> >();
    test_atom_expr<"[r2]",  RingSize<2> >();
    test_atom_expr<"[r42]", RingSize<42> >();
}

TEST(TestCTSmarts, RingConnectivity)
{
    test_atom_expr<"[x1]",  RingConnectivity<1> >();
    test_atom_expr<"[x2]",  RingConnectivity<2> >();
    test_atom_expr<"[x42]", RingConnectivity<42> >();
}


TEST(TestCTSmarts, Charge)
{
    test_atom_expr<"[-]",  Charge<-1> >();
    test_atom_expr<"[-0]", Charge<0> >();
    test_atom_expr<"[-1]", Charge<-1> >();
    test_atom_expr<"[-3]", Charge<-3> >();

    test_atom_expr<"[+]",  Charge<1> >();
    test_atom_expr<"[+0]", Charge<0> >();
    test_atom_expr<"[+1]", Charge<1> >();
    test_atom_expr<"[+3]", Charge<3> >();

    test_atom_expr<"[++]",  Charge<2> >();
    test_atom_expr<"[+++]",  Charge<3> >();
    test_atom_expr<"[--]",  Charge<-2> >();
    test_atom_expr<"[---]",  Charge<-3> >();
}

TEST(TestCTSmarts, Chiral)
{
}

TEST(TestCTSmarts, Class)
{
}

TEST(TestCTSmarts, BracketAtom)
{
    using N = AliphaticAtom<7>;
    using O = AliphaticAtom<8>;
    using F = AliphaticAtom<9>;

    test_atom_expr<"[N]", N >();

    // operations
    test_atom_expr<"[FN]",    And<F, N> >();
    test_atom_expr<"[FN,O]",  Or<And<F, N>, O> >();
    test_atom_expr<"[!N]",    Not<N> >();
    test_atom_expr<"[!!N]",   N >();
    test_atom_expr<"[!!!N]",  Not<N> >();
    test_atom_expr<"[!!!!N]", N >();
    test_atom_expr<"[O!N]",   And<O, Not<N>> >();
    test_atom_expr<"[F,N,O]",  Or<F, N, O> >();

    // isotope


    // element

}


TEST(TestCTSmarts, BondPrimitive)
{
    test_bond_expr<"**",    ImplicitBond >();
    test_bond_expr<"*-*",   BondOrder<1> >();
    test_bond_expr<"*=*",   BondOrder<2> >();
    test_bond_expr<"*#*",   BondOrder<3> >();
    test_bond_expr<"*$*",   BondOrder<4> >();
    test_bond_expr<"*:*",   AromaticBond >();
    test_bond_expr<"*~*",   AnyBond >();
    test_bond_expr<"*@*",   RingBond >();

    //test_bond_expr<"*/*",   UpBond >();
    //test_bond_expr<"*\\*",  DownBond >();
    //test_bond_expr<"*/?*",  UpOrDownBond >();
    //test_bond_expr<"*\\?*", UpOrDownBond >();


    test_bond_expr<"[*][*]",  ImplicitBond >();
    test_bond_expr<"[*]-[*]",  BondOrder<1> >();
    test_bond_expr<"[*]=[*]",  BondOrder<2> >();
}

TEST(TestCTSmarts, BondExpr)
{
    using Single = BondOrder<1>;
    using Double = BondOrder<2>;

    test_bond_expr<"*!-*", Not<Single> >();
    test_bond_expr<"*-,=*", Or<Single, Double> >();
    test_bond_expr<"*-@*", And<Single, RingBond> >();
    test_bond_expr<"*-&@*", And<Single, RingBond> >();
    test_bond_expr<"*-;@*", And<Single, RingBond> >();
    test_bond_expr<"*-@,=@*", Or<And<Single, RingBond>, And<Double, RingBond>> >();
    test_bond_expr<"*-&@,=&@*", Or<And<Single, RingBond>, And<Double, RingBond>> >();
    test_bond_expr<"*-,=;@*", And<Or<Single, Double>, RingBond> >();
}

TEST(TestCTSmarts, RingBond)
{
    using Single = BondOrder<1>;
    using Double = BondOrder<2>;

    test_bond<"*1***1", Bond<3, 3, 0, ImplicitBond>, 3 >();
    test_bond<"*-1***1", Bond<3, 3, 0, Single>, 3 >();
    test_bond<"*1***-1", Bond<3, 3, 0, Single>, 3 >();
    test_bond<"*=1***=1", Bond<3, 3, 0, Double>, 3 >();


    //         0  2   4
    test_bond<"*1**12**2", Bond<2, 2, 0, ImplicitBond>, 2 >();
    test_bond<"*1**12**2", Bond<5, 4, 2, ImplicitBond>, 5 >();
    test_bond<"*2**21**1", Bond<2, 2, 0, ImplicitBond>, 2 >();
    test_bond<"*2**21**1", Bond<5, 4, 2, ImplicitBond>, 5 >();
    test_bond<"*1**11**1", Bond<2, 2, 0, ImplicitBond>, 2 >();
    test_bond<"*1**11**1", Bond<5, 4, 2, ImplicitBond>, 5 >();

    // adjacent ring bonds
    test_bond<"*:1=2**:1*=2", Bond<2, 2, 0, AromaticBond>, 2 >();
    test_bond<"*:1=2**:1*=2", Bond<4, 3, 0, Double>, 4 >();
    test_bond<"*:1=2**1*2", Bond<2, 2, 0, AromaticBond>, 2 >();
    test_bond<"*:1=2**1*2", Bond<4, 3, 0, Double>, 4 >();
    test_bond<"*12**:1*=2", Bond<2, 2, 0, AromaticBond>, 2 >();
    test_bond<"*12**:1*=2", Bond<4, 3, 0, Double>, 4 >();


    // conflicing ring bonds
    test_error<"*-1***=1", ConflicingRingBondError>();
}

TEST(TestCTSmarts, Operators)
{
    using C = AliphaticAtom<6>;
    using N = AliphaticAtom<7>;
    using O = AliphaticAtom<8>;
    using F = AliphaticAtom<9>;


    test_atom_expr<"[C&N]", And<C, N> >();
    test_atom_expr<"[C,N]",      Or<C, N> >();
    test_atom_expr<"[C;N]",  And<C, N> >();

    test_atom_expr<"[!CN]",  And<Not<C>, N> >();
    test_atom_expr<"[!C&N]", And<Not<C>, N> >();
    test_atom_expr<"[!C,N]",      Or<Not<C>, N> >();
    test_atom_expr<"[!C;N]",  And<Not<C>, N> >();

    test_atom_expr<"[C!N]",  And<C, Not<N>> >();
    test_atom_expr<"[C&!N]", And<C, Not<N>> >();
    test_atom_expr<"[C,!N]",      Or<C, Not<N>> >();
    test_atom_expr<"[C;!N]",  And<C, Not<N>> >();

    test_atom_expr<"[!C!N]",  And<Not<C>, Not<N>> >();
    test_atom_expr<"[!C&!N]", And<Not<C>, Not<N>> >();
    test_atom_expr<"[!C,!N]",      Or<Not<C>, Not<N>> >();
    test_atom_expr<"[!C;!N]",  And<Not<C>, Not<N>> >();



    test_atom_expr<"[CN]",    And<C, N> >();
    test_atom_expr<"[CNO]",   And<C, N, O> >();
    test_atom_expr<"[C,NO]",  Or<C, And<N, O>> >();
    test_atom_expr<"[C;NO]",  And<C, And<N, O>> >();

    test_atom_expr<"[C,N]",    Or<C, N> >();
    test_atom_expr<"[CN,O]",   Or<And<C, N>, O> >();
    test_atom_expr<"[C,N,O]",  Or<C, N, O> >();
    test_atom_expr<"[C;N,O]",  And<C, Or<N, O>> >();

    test_atom_expr<"[C;N]",    And<C, N> >();
    test_atom_expr<"[CN;O]",   And<And<C, N>, O> >();
    test_atom_expr<"[C,N;O]",  And<Or<C, N>, O> >();
    test_atom_expr<"[C;N;O]",  And<C, N, O> >();




    test_atom_expr<"[CNOF]",  And<C, N, O, F> >();

    test_atom_expr<"[C,NOF]",  Or<C, And<N, O, F>> >();
    test_atom_expr<"[CN,OF]",  Or<And<C, N>, And<O, F>> >();
    test_atom_expr<"[CNO,F]",  Or<And<C, N, O>, F> >();

    test_atom_expr<"[C;NOF]",  And<C, And<N, O, F>> >();
    test_atom_expr<"[CN;OF]",  And<And<C, N>, And<O, F>> >();
    test_atom_expr<"[CNO;F]",  And<And<C, N, O>, F> >();

    test_atom_expr<"[C,N;OF]",  And<Or<C, N>, And<O, F>> >();
    test_atom_expr<"[C,NO;F]",  And<Or<C, And<N, O>>, F> >();
    test_atom_expr<"[C;N,OF]",  And<C, Or<N, And<O, F>>> >();
    test_atom_expr<"[C;NO,F]",  And<C, Or<And<N, O>, F>> >();
    test_atom_expr<"[CN,O;F]",  And<Or<And<C, N>, O>, F> >();
    test_atom_expr<"[CN;O,F]",  And<And<C, N>, Or<O, F>> >();




    /*
    test_atom_expr<"[CN]**",  AndHigh<C, N> >();
    test_atom_expr<"[C&N]**", AndHigh<C, N> >();
    test_atom_expr<"[C,N]**", Or<C, N> >();
    test_atom_expr<"[C;N]**", And<C, N> >();
    */
}

TEST(TestCTSmarts, GetCentralAtom)
{
    EXPECT_EQ(getCentralAtom(Smarts<"C">{}), -1);
    EXPECT_EQ(getCentralAtom(Smarts<"CC">{}), -1);
    EXPECT_EQ(getCentralAtom(Smarts<"CCC">{}), 1);
    EXPECT_EQ(getCentralAtom(Smarts<"C(C)C">{}), 0);
    EXPECT_EQ(getCentralAtom(Smarts<"CCCC">{}), -1);
    EXPECT_EQ(getCentralAtom(Smarts<"CC(C)C">{}), 1);
    EXPECT_EQ(getCentralAtom(Smarts<"C(C)(C)C">{}), 0);
    EXPECT_EQ(getCentralAtom(Smarts<"CCCCC">{}), -1);
    EXPECT_EQ(getCentralAtom(Smarts<"CC(C)CC">{}), -1);
    EXPECT_EQ(getCentralAtom(Smarts<"CC(C)(C)C">{}), 1);
    EXPECT_EQ(getCentralAtom(Smarts<"C(C)(C)(C)C">{}), 0);
    EXPECT_EQ(getCentralAtom(Smarts<"C1CC1">{}), -1);
}


TEST(TestCTSmarts, MACCS)
{
    /*
    static_assert(Smarts<"[#103,#104,#105,#106,#107,#106,#109,#110,#111,#112]">::valid);
    static_assert(Smarts<"[#103,#104]">::valid);
    static_assert(Smarts<"[Ge,As,Se,Sn,Sb,Te,Tl,Pb,Bi]">::valid);
    static_assert(Smarts<"[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]">::valid);
    static_assert(Smarts<"[Sc,Ti,Y,Zr,Hf]">::valid);
    static_assert(Smarts<"[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]">::valid);
    static_assert(Smarts<"[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]">::valid);
    static_assert(Smarts<"[!#6;!#1]1~*~*~*~1">::valid);
    static_assert(Smarts<"[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]">::valid);
    static_assert(Smarts<"[Be,Mg,Ca,Sr,Ba,Ra]">::valid);

    static_assert(Smarts<"*1~*~*~*~1">::valid);
    static_assert(Smarts<"[Cu,Zn,Ag,Cd,Au,Hg]">::valid);
    static_assert(Smarts<"[#8]~[#7](~[#6])~[#6]">::valid);
    static_assert(Smarts<"[#16]-[#16]">::valid);
    static_assert(Smarts<"[#8]~[#6](~[#8])~[#8]">::valid);
    static_assert(Smarts<"[!#6;!#1]1~*~*~1">::valid);
    static_assert(Smarts<"[#6]#[#6]">::valid);
    static_assert(Smarts<"[B,Al,Ga,In,Tl]">::valid);
    static_assert(Smarts<"*1~*~*~*~*~*~*~1">::valid);
    static_assert(Smarts<"[Si]">::valid);


    static_assert(Smarts<"*!:*:*!:*">::valid);
    static_assert(Smarts<"*!@*@*!@*">::valid);
    static_assert(Smarts<"*!@[#7]@*">::valid);
    static_assert(Smarts<"*!@[#8]!@*">::valid);
    static_assert(Smarts<"*!@[CH2]!@*">::valid);
    static_assert(Smarts<"*1~*~*~*~*~*~1">::valid);
    static_assert(Smarts<"*1~*~*~*~*~*~1">::valid);
    static_assert(Smarts<"*1~*~*~*~*~1">::valid);
    static_assert(Smarts<"*1~*~*~1">::valid);
    static_assert(Smarts<"*@*!@*@*">::valid);
    static_assert(Smarts<"*@*!@[#16]">::valid);
    static_assert(Smarts<"*@*!@[#7]">::valid);
    static_assert(Smarts<"*@*!@[#8]">::valid);
    static_assert(Smarts<"*@*!@[#8]">::valid);
    static_assert(Smarts<"*@*(@*)@*">::valid);
    static_assert(Smarts<"*~*(~*)(~*)~*">::valid);
    static_assert(Smarts<"*~[!#6;!#1](~*)~*">::valid);
    static_assert(Smarts<"*~[#16](~*)~*">::valid);
    static_assert(Smarts<"*~[#7](~*)~*">::valid);
    static_assert(Smarts<"*~[CH2]~[!#6;!#1;!H0]">::valid);
    static_assert(Smarts<"*~[CH2]~[#7]">::valid);
    static_assert(Smarts<"*~[CH2]~[#8]">::valid);
    static_assert(Smarts<"Br">::valid);
    static_assert(Smarts<"Cl">::valid);
    static_assert(Smarts<"F">::valid);
    static_assert(Smarts<"[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]">::valid);
    static_assert(Smarts<"[!#6;!#1;!H0]">::valid);
    static_assert(Smarts<"[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]">::valid);
    static_assert(Smarts<"[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]">::valid);
    static_assert(Smarts<"[!#6;!#1;!H0]~*~[CH2]~*">::valid);
    static_assert(Smarts<"[!#6;!#1;!H0]~[!#6;!#1;!H0]">::valid);
    static_assert(Smarts<"[!#6;!#1]1~*~*~*~*~*~1">::valid);
    static_assert(Smarts<"[!#6;!#1]1~*~*~*~*~1">::valid);
    static_assert(Smarts<"[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[!#6;!#1;!H0]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[!#6;!#1]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[!#6;!#1]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[#16]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[#16]~[!#6;!#1]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[#7]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[#7]~[!#6;!#1]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[#8]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[CH2]~*">::valid);
    static_assert(Smarts<"[!#6;!#1]~[CH2]~*">::valid);
    static_assert(Smarts<"[!#6;!#1]~[CH2]~[!#6;!#1]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[CH3]">::valid);
    static_assert(Smarts<"[!#6;!#1]~[F,Cl,Br,I]">::valid);
    static_assert(Smarts<"[!#6;R]">::valid);
    //static_assert(Smarts<"[!+0]">::valid);
    static_assert(Smarts<"[!C;!c;!#1;!H0]~*~[!C;!c;!#1;!H0]">::valid);
    static_assert(Smarts<"[!C;!c;R]">::valid);
    static_assert(Smarts<"[#15]">::valid);
    static_assert(Smarts<"[#16R]">::valid);
    static_assert(Smarts<"[#16]">::valid);
    static_assert(Smarts<"[#16]!:*:*">::valid);
    static_assert(Smarts<"[#16]-[#8]">::valid);
    static_assert(Smarts<"[#16]=*">::valid);
    static_assert(Smarts<"[#16]=[#8]">::valid);
    static_assert(Smarts<"[#16]~*(~*)~*">::valid);
    static_assert(Smarts<"[#16]~*~[#7]">::valid);
    static_assert(Smarts<"[#6]#[#7]">::valid);
    static_assert(Smarts<"[#6]-[#7]">::valid);
    static_assert(Smarts<"[#6]-[#8]">::valid);
    static_assert(Smarts<"[#6]=;@[#6](@*)@*">::valid);
    static_assert(Smarts<"[#6]=[#6]">::valid);
    static_assert(Smarts<"[#6]=[#6](~*)~*">::valid);
    static_assert(Smarts<"[#6]=[#6](~[!#6;!#1])~[!#6;!#1]">::valid);
    static_assert(Smarts<"[#6]=[#6](~[#6])~[#6]">::valid);
    static_assert(Smarts<"[#6]=[#6]~[#7]">::valid);
    static_assert(Smarts<"[#6]=[#7]">::valid);
    static_assert(Smarts<"[#6]=[#8]">::valid);
    static_assert(Smarts<"[#6]~[!#6;!#1](~[#6])(~[#6])~*">::valid);
    static_assert(Smarts<"[#6]~[#16]~[#7]">::valid);
    static_assert(Smarts<"[#6]~[#16]~[#8]">::valid);
    static_assert(Smarts<"[#6]~[#6](~[#6])(~[#6])~*">::valid);
    static_assert(Smarts<"[#6]~[#7](~[#6])~[#6]">::valid);
    static_assert(Smarts<"[#7;!H0]">::valid);
    static_assert(Smarts<"[#7;R]">::valid);
    static_assert(Smarts<"[#7]">::valid);
    static_assert(Smarts<"[#7]!:*:*">::valid);
    static_assert(Smarts<"[#7]">::valid);
    static_assert(Smarts<"[#7]-[#8]">::valid);
    static_assert(Smarts<"[#7]=*">::valid);
    static_assert(Smarts<"[#7]=[#8]">::valid);
    static_assert(Smarts<"[#7]~*(~*)~*">::valid);
    static_assert(Smarts<"[#7]~*~*~*~[#7]">::valid);
    static_assert(Smarts<"[#7]~*~*~*~[#8]">::valid);
    static_assert(Smarts<"[#7]~*~*~[#7]">::valid);
    static_assert(Smarts<"[#7]~*~*~[#8]">::valid);
    static_assert(Smarts<"[#7]~*~[#7]">::valid);
    static_assert(Smarts<"[#7]~*~[#8]">::valid);
    static_assert(Smarts<"[#7]~*~[CH2]~*">::valid);
    static_assert(Smarts<"[#7]~[#16]">::valid);
    static_assert(Smarts<"[#7]~[#6](~[#6])~[#7]">::valid);
    static_assert(Smarts<"[#7]~[#6](~[#7])~[#7]">::valid);
    static_assert(Smarts<"[#7]~[#6](~[#8])~[#7]">::valid);
    static_assert(Smarts<"[#7]~[#6](~[#8])~[#8]">::valid);
    static_assert(Smarts<"[#7]~[#6]~[#8]">::valid);
    static_assert(Smarts<"[#7]~[#7]">::valid);
    static_assert(Smarts<"[#7]~[#8]">::valid);
    static_assert(Smarts<"[#8R]">::valid);
    static_assert(Smarts<"[#8]">::valid);
    static_assert(Smarts<"[#8]!:*:*">::valid);
    static_assert(Smarts<"[#8]">::valid);
    static_assert(Smarts<"[#8]">::valid);
    static_assert(Smarts<"[#8]">::valid);
    static_assert(Smarts<"[#8]=*">::valid);
    static_assert(Smarts<"[#8]~*~*~*~[#8]">::valid);
    static_assert(Smarts<"[#8]~*~*~[#8]">::valid);
    static_assert(Smarts<"[#8]~*~[CH2]~*">::valid);
    static_assert(Smarts<"[#8]~[!#6;!#1](~[#8])(~[#8])">::valid);
    static_assert(Smarts<"[#8]~[#16](~[#8])~[#8]">::valid);
    static_assert(Smarts<"[#8]~[#16]~[#8]">::valid);
    static_assert(Smarts<"[#8]~[#6](~[#6])~[#6]">::valid);
    static_assert(Smarts<"[#8]~[#6](~[#7])~[#6]">::valid);
    static_assert(Smarts<"[#8]~[#6]~[#8]">::valid);
    static_assert(Smarts<"[#8]~[#7](~[#8])~[#6]">::valid);
//    static_assert(Smarts<"[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]">::valid);
//    static_assert(Smarts<"[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]">::valid);
//    static_assert(Smarts<"[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]">::valid);
//    static_assert(Smarts<"[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]">::valid);
//    static_assert(Smarts<"[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]">::valid);
//    static_assert(Smarts<"[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]">::valid);
//    static_assert(Smarts<"[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]">::valid);
//    static_assert(Smarts<"[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]">::valid);
    static_assert(Smarts<"[C;H2,H3][!#6;!#1][C;H2,H3]">::valid);
    static_assert(Smarts<"[C;H3,H4]">::valid);
    static_assert(Smarts<"[C;H3,H4]">::valid);
    static_assert(Smarts<"[CH2]=*">::valid);
    static_assert(Smarts<"[CH3]">::valid);
    static_assert(Smarts<"[CH3]~*~*~*~[CH2]~*">::valid);
    static_assert(Smarts<"[CH3]~*~[CH2]~*">::valid);
    static_assert(Smarts<"[CH3]~*~[CH3]">::valid);
    static_assert(Smarts<"[CH3]~[CH2]~*">::valid);
    static_assert(Smarts<"[F,Cl,Br,I]">::valid);
    static_assert(Smarts<"[F,Cl,Br,I]!@*@*">::valid);
    static_assert(Smarts<"[F,Cl,Br,I]~*(~*)~*">::valid);
    static_assert(Smarts<"[I]">::valid);
    static_assert(Smarts<"[Li,Na,K,Rb,Cs,Fr]">::valid);
    static_assert(Smarts<"[NH2]">::valid);
    static_assert(Smarts<"[O;!H0]">::valid);
    static_assert(Smarts<"[R]">::valid);
    static_assert(Smarts<"a">::valid);
    static_assert(Smarts<"c:n">::valid);
    */
}

auto mockAcetateAnion()
{
    // SMILES: CC(=O)[O-]
    // Generated by SmiToMock
    auto mol = Kitimar::Molecule::MockMolecule{
        {{6,0,0,1,3,false,false},{6,0,0,3,0,false,false},{8,0,0,1,0,false,false},{8,0,-1,1,0,false,false}},
        {{0,1,1,false,false},{1,2,2,false,false},{1,3,1,false,false}}
    };
    return mol;
}

auto mockButane()
{
    // SMILES: CCCC
    // Generated by SmiToMock
    auto mol = Kitimar::Molecule::MockMolecule{
        {{6,0,0,1,3,false,false},{6,0,0,2,2,false,false},{6,0,0,2,2,false,false},{6,0,0,1,3,false,false}},
        {{0,1,1,false,false},{1,2,1,false,false},{2,3,1,false,false}}
    };
    return mol;
}



TEST(TestCTSmarts, CTSmarts_contains)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]

    // single atom
    EXPECT_TRUE(CTSmarts::contains<"C">(mol));
    EXPECT_TRUE(CTSmarts::contains<"O">(mol));
    EXPECT_TRUE(CTSmarts::contains<"[O-]">(mol));

    EXPECT_FALSE(CTSmarts::contains<"N">(mol));
    EXPECT_FALSE(CTSmarts::contains<"[O+]">(mol));

    // single bond
    EXPECT_TRUE(CTSmarts::contains<"CC">(mol));
    EXPECT_TRUE(CTSmarts::contains<"C=O">(mol));
    EXPECT_TRUE(CTSmarts::contains<"C[O-]">(mol));

    EXPECT_FALSE(CTSmarts::contains<"C=C">(mol));
    EXPECT_FALSE(CTSmarts::contains<"C#O">(mol));
    EXPECT_FALSE(CTSmarts::contains<"C[O+]">(mol));

    // general case
    EXPECT_TRUE(CTSmarts::contains<"CC(=O)[O-]">(mol));
    EXPECT_TRUE(CTSmarts::contains<"CC([O-])=O">(mol));
    EXPECT_TRUE(CTSmarts::contains<"O=C(C)[O-]">(mol));

    EXPECT_FALSE(CTSmarts::contains<"CC(=O)N">(mol));
}

TEST(TestCTSmarts, CTSmarts_atom)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]

    auto C0 = get_atom(mol, 0);
    auto C1 = get_atom(mol, 1);
    auto O2 = get_atom(mol, 2);
    auto O3 = get_atom(mol, 3);

    // single atom
    EXPECT_TRUE(CTSmarts::atom<"C">(mol, C0));
    EXPECT_TRUE(CTSmarts::atom<"C">(mol, C1));
    EXPECT_TRUE(CTSmarts::atom<"O">(mol, O2));
    EXPECT_TRUE(CTSmarts::atom<"[O-]">(mol, O3));

    EXPECT_FALSE(CTSmarts::atom<"C">(mol, O2));
    EXPECT_FALSE(CTSmarts::atom<"[O-]">(mol, O2));

    // single bond
    EXPECT_TRUE(CTSmarts::atom<"CC">(mol, C0));
    EXPECT_TRUE(CTSmarts::atom<"CC">(mol, C1));
    EXPECT_TRUE(CTSmarts::atom<"CO">(mol, C1));
    EXPECT_TRUE(CTSmarts::atom<"O=C">(mol, O2));
    EXPECT_TRUE(CTSmarts::atom<"OC">(mol, O3));

    EXPECT_FALSE(CTSmarts::atom<"CO">(mol, C0));
    EXPECT_FALSE(CTSmarts::atom<"CO">(mol, O2));
    EXPECT_FALSE(CTSmarts::atom<"OC">(mol, C1));

    // general case
    EXPECT_TRUE(CTSmarts::atom<"CC(=O)[O-]">(mol, C0));
    EXPECT_TRUE(CTSmarts::atom<"C(C)(=O)[O-]">(mol, C1));
    EXPECT_TRUE(CTSmarts::atom<"O=CO">(mol, O2));
    EXPECT_TRUE(CTSmarts::atom<"OC=O">(mol, O3));

    EXPECT_FALSE(CTSmarts::atom<"CC(=O)[O-]">(mol, C1));
    EXPECT_FALSE(CTSmarts::atom<"C(C)(=O)[O-]">(mol, C0));
    EXPECT_FALSE(CTSmarts::atom<"O=CO">(mol, O3));
    EXPECT_FALSE(CTSmarts::atom<"OC=O">(mol, O2));
}

TEST(TestCTSmarts, CTSmarts_bond)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]

    auto CC0 = get_bond(mol, 0);
    auto CO1 = get_bond(mol, 1);
    auto CO2 = get_bond(mol, 2);

    // single bond
    EXPECT_TRUE(CTSmarts::bond<"CC">(mol, CC0));
    EXPECT_TRUE(CTSmarts::bond<"C=O">(mol, CO1));
    EXPECT_TRUE(CTSmarts::bond<"O=C">(mol, CO1));
    EXPECT_TRUE(CTSmarts::bond<"CO">(mol, CO2));
    EXPECT_TRUE(CTSmarts::bond<"OC">(mol, CO2));

    EXPECT_FALSE(CTSmarts::bond<"C=C">(mol, CC0));
    EXPECT_FALSE(CTSmarts::bond<"NC">(mol, CC0));
    EXPECT_FALSE(CTSmarts::bond<"CN">(mol, CC0));
    EXPECT_FALSE(CTSmarts::bond<"CO">(mol, CO1));
    EXPECT_FALSE(CTSmarts::bond<"N=O">(mol, CO1));
    EXPECT_FALSE(CTSmarts::bond<"C=N">(mol, CO1));
    EXPECT_FALSE(CTSmarts::bond<"C=O">(mol, CO2));

    // general case
    EXPECT_TRUE(CTSmarts::bond<"CCO">(mol, CC0));
    EXPECT_TRUE(CTSmarts::bond<"C(C)O">(mol, CC0));
    EXPECT_TRUE(CTSmarts::bond<"O=CO">(mol, CO1));
    EXPECT_TRUE(CTSmarts::bond<"C(=O)O">(mol, CO1));
    EXPECT_TRUE(CTSmarts::bond<"OC=O">(mol, CO2));
    EXPECT_TRUE(CTSmarts::bond<"C(O)=O">(mol, CO2));

    EXPECT_FALSE(CTSmarts::bond<"CCO">(mol, CO1));
    EXPECT_FALSE(CTSmarts::bond<"OCC">(mol, CC0));
    EXPECT_FALSE(CTSmarts::bond<"O(C)C">(mol, CC0));
    EXPECT_FALSE(CTSmarts::bond<"OC=O">(mol, CO1));
    EXPECT_FALSE(CTSmarts::bond<"C(O)=O">(mol, CO1));
    EXPECT_FALSE(CTSmarts::bond<"O=CO">(mol, CO2));
    EXPECT_FALSE(CTSmarts::bond<"C(=O)O">(mol, CO2));
}


void compare_maps(auto &&maps, const auto &refMaps)
{

    std::cout << "maps:" << std::endl;
    for (const auto &map : maps)
        std::cout << "    " << map << std::endl;

    std::cout << "refMaps:" << std::endl;
    for (const auto &map : refMaps)
        std::cout << "    " << map << std::endl;


    auto map = std::begin(maps);
    auto refMap = std::begin(refMaps);
    while (map != std::end(maps) && refMap != std::end(refMaps)) {
        EXPECT_TRUE(std::ranges::equal(*map, *refMap));
        ++map;
        ++refMap;
    }
    EXPECT_EQ(map, std::end(maps));
    EXPECT_EQ(refMap, std::end(refMaps));
}

template<ctll::fixed_string SMARTS>
void test_unique(auto &mol, std::initializer_list<std::array<int, Smarts<SMARTS>::numAtoms>> refMaps)
{
    compare_maps(CTSmarts::multi<SMARTS>(mol), refMaps);
}

template<ctll::fixed_string SMARTS>
void test_all(auto &mol, std::initializer_list<std::array<int, Smarts<SMARTS>::numAtoms>> refMaps)
{
    compare_maps(CTSmarts::multi<SMARTS>(mol, CTSmarts::All), refMaps);
}

TEST(TestCTSmarts, CTSmarts_multi)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]
    auto mol2 = mockButane(); // CCCCC

    test_unique<"*~*">(mol, { {0, 1}, {1, 2}, {1, 3} });
    test_unique<"*~*~*">(mol, { {0, 1, 2}, {0, 1, 3}, {2, 1, 3} });
    test_unique<"*~*~*~*">(mol2, { {0, 1, 2, 3} });
    test_unique<"*~*(~*)~*">(mol, { {0, 1, 2, 3} });

    test_all<"*~*">(mol, { {0, 1}, {1, 0}, {1, 2}, {2, 1}, {1, 3}, {3, 1} });
    test_all<"*~*~*">(mol, { {0, 1, 2}, {0, 1, 3}, {2, 1, 0}, {2, 1, 3}, {3, 1, 0}, {3, 1, 2} });
    test_all<"*~*~*~*">(mol2, { {0, 1, 2, 3}, {3, 2, 1, 0} });
    test_all<"*~*(~*)~*">(mol, { {0, 1, 2, 3}, {0, 1, 3, 2}, {2, 1, 0, 3}, {2, 1, 3, 0}, {3, 1, 0, 2}, {3, 1, 2, 0} });
}

TEST(TestCTSmarts, CTSmarts_capture)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]

    auto C0 = get_atom(mol, 0);
    auto C1 = get_atom(mol, 1);
    auto O2 = get_atom(mol, 2);
    auto O3 = get_atom(mol, 3);

    // single atom
    EXPECT_EQ(CTSmarts::capture<"C">(mol), std::make_tuple(true, C0));
    EXPECT_EQ(CTSmarts::capture<"O">(mol), std::make_tuple(true, O2));

    EXPECT_EQ(CTSmarts::capture<"N">(mol), std::make_tuple(false, null_atom(mol)));

    // single bond
    EXPECT_EQ(CTSmarts::capture<"CC">(mol), std::make_tuple(true, C0, C1));
    EXPECT_EQ(CTSmarts::capture<"C=O">(mol), std::make_tuple(true, C1, O2));
    EXPECT_EQ(CTSmarts::capture<"CO">(mol), std::make_tuple(true, C1, O3));

    EXPECT_EQ(CTSmarts::capture<"[C:1]=[O:2]">(mol), std::make_tuple(true, C1, O2));
    EXPECT_EQ(CTSmarts::capture<"[C:2]=[O:1]">(mol), std::make_tuple(true, O2, C1));
    EXPECT_EQ(CTSmarts::capture<"[C:1]=O">(mol), std::make_tuple(true, C1));
    EXPECT_EQ(CTSmarts::capture<"C=[O:1]">(mol), std::make_tuple(true, O2));

    EXPECT_EQ(CTSmarts::capture<"*#*">(mol), std::make_tuple(false, null_atom(mol), null_atom(mol)));

    // general case
    EXPECT_EQ(CTSmarts::capture<"CC=O">(mol), std::make_tuple(true, C0, C1, O2));
    EXPECT_EQ(CTSmarts::capture<"O=CC">(mol), std::make_tuple(true, O2, C1, C0));
    EXPECT_EQ(CTSmarts::capture<"C(=O)C">(mol), std::make_tuple(true, C1, O2, C0));
    EXPECT_EQ(CTSmarts::capture<"CC(=O)O">(mol), std::make_tuple(true, C0, C1, O2, O3));
    EXPECT_EQ(CTSmarts::capture<"CC(O)=O">(mol), std::make_tuple(true, C0, C1, O3, O2));
}

TEST(TestCTSmarts, CTSmarts_captureAtom)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]

    auto C0 = get_atom(mol, 0);
    auto C1 = get_atom(mol, 1);
    auto O2 = get_atom(mol, 2);
    auto O3 = get_atom(mol, 3);

    // single atom
    EXPECT_EQ(CTSmarts::captureAtom<"C">(mol), C0);
    EXPECT_EQ(CTSmarts::captureAtom<"O">(mol), O2);

    EXPECT_EQ(CTSmarts::captureAtom<"N">(mol), null_atom(mol));

    // single bond
    EXPECT_EQ(CTSmarts::captureAtom<"CC">(mol), C0);
    EXPECT_EQ(CTSmarts::captureAtom<"C=O">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"CO">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"O=C">(mol), O2);
    EXPECT_EQ(CTSmarts::captureAtom<"OC">(mol), O3);

    EXPECT_EQ(CTSmarts::captureAtom<"[C:1]C">(mol), C0);
    EXPECT_EQ(CTSmarts::captureAtom<"[C:1]=O">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"[C:1]O">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"[O:1]=C">(mol), O2);
    EXPECT_EQ(CTSmarts::captureAtom<"[O:1]C">(mol), O3);

    EXPECT_EQ(CTSmarts::captureAtom<"C[C:1]">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"C=[O:1]">(mol), O2);
    EXPECT_EQ(CTSmarts::captureAtom<"C[O:1]">(mol), O3);
    EXPECT_EQ(CTSmarts::captureAtom<"O=[C:1]">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"O[C:1]">(mol), C1);

    EXPECT_EQ(CTSmarts::captureAtom<"C#C">(mol), null_atom(mol));

    // general case
    EXPECT_EQ(CTSmarts::captureAtom<"CCO">(mol), C0);
    EXPECT_EQ(CTSmarts::captureAtom<"CO">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"OCC">(mol), O3);
    EXPECT_EQ(CTSmarts::captureAtom<"O=CC">(mol), O2);
    EXPECT_EQ(CTSmarts::captureAtom<"CC(=O)O">(mol), C0);
    EXPECT_EQ(CTSmarts::captureAtom<"C(C)(=O)O">(mol), C1);

    EXPECT_EQ(CTSmarts::captureAtom<"[C:1]C(=O)O">(mol), C0);
    EXPECT_EQ(CTSmarts::captureAtom<"C[C:1](=O)O">(mol), C1);
    EXPECT_EQ(CTSmarts::captureAtom<"CC(=[O:1])O">(mol), O2);
    EXPECT_EQ(CTSmarts::captureAtom<"CC(=O)[O:1]">(mol), O3);

    EXPECT_EQ(CTSmarts::captureAtom<"CC(=O)N">(mol), null_atom(mol));
}

template<typename T>
struct identify_type;

TEST(TestCTSmarts, EdgeList)
{
    auto smarts = Smarts<"*1*(*)*1">{};

    constexpr auto edges = EdgeList{smarts}.data;

    static_assert(edges.size() == 4);
    static_assert(edges[0] == Edge{0, 0, 1});
    static_assert(edges[1] == Edge{1, 1, 2});
    static_assert(edges[2] == Edge{2, 1, 3});
    static_assert(edges[3] == Edge{3, 3, 0});
}


TEST(TestCTSmarts, DfsSearch)
{

    //auto smarts = Smarts<"CC(C)C">{};
    //auto smarts = Smarts<"*1**1">{};
    //auto smarts = Smarts<"*1**1*">{};
    auto smarts = Smarts<"*1**12**2*">{};

    auto edges = EdgeList(smarts);
    auto degrees = VertexDegree(smarts, edges);
    auto adjList = IncidentList(smarts, edges, degrees);

    //DfsSearchEventsVisitor visitor(smarts);
    //dfsSearch(smarts, visitor, adjList);

    auto dfsSearchEvents = DfsSearchEvents(smarts, adjList);

    //for (const auto &event : visitor.events)
    for (const auto &event : dfsSearchEvents.events)
        std::cout << event << std::endl;

}


TEST(TestCTSmarts, ExpressionFrequency)
{
    using C = AliphaticAtom<6>;
    using O = AliphaticAtom<8>;
    //using S = AliphaticAtom<16>;


    static_assert(expressionFrequency(Not<C>{}) == 1 - expressionFrequency(C{}));

    static_assert(expressionFrequency(And<C, O>{}) == expressionFrequency(O{}));
    static_assert(expressionFrequency(And<O, C>{}) == expressionFrequency(O{}));

    static_assert(expressionFrequency(Or<C, O>{}) == expressionFrequency(C{}));
    static_assert(expressionFrequency(Or<O, C>{}) == expressionFrequency(C{}));




}



template<typename Compare, typename Expr, typename OptimizedExpr>
constexpr void test_selection_sort(Expr expr, OptimizedExpr)
{
    static_assert(std::is_same_v<decltype(selectionSort<ProjExprFrequency, Compare>(expr)), OptimizedExpr>);
}

TEST(TestCTSmarts, CtllSort)
{
    using C = AliphaticAtom<6>;
    using O = AliphaticAtom<8>;
    using S = AliphaticAtom<16>;

    auto less = ctll::list<S, O, C>{};
    auto greater = ctll::list<C, O, S>{};

    auto permutation1 = ctll::list<C, O, S>{};
    auto permutation2 = ctll::list<C, S, O>{};
    auto permutation3 = ctll::list<O, C, S>{};
    auto permutation4 = ctll::list<O, S, C>{};
    auto permutation5 = ctll::list<S, C, O>{};
    auto permutation6 = ctll::list<S, O, C>{};

    test_selection_sort<std::less<>>(permutation1, less);
    test_selection_sort<std::less<>>(permutation2, less);
    test_selection_sort<std::less<>>(permutation3, less);
    test_selection_sort<std::less<>>(permutation4, less);
    test_selection_sort<std::less<>>(permutation5, less);
    test_selection_sort<std::less<>>(permutation6, less);

    test_selection_sort<std::greater<>>(permutation1, greater);
    test_selection_sort<std::greater<>>(permutation2, greater);
    test_selection_sort<std::greater<>>(permutation3, greater);
    test_selection_sort<std::greater<>>(permutation4, greater);
    test_selection_sort<std::greater<>>(permutation5, greater);
    test_selection_sort<std::greater<>>(permutation6, greater);
}

/*
template<typename Compare, typename Expr, typename OptimizedExpr>
constexpr void test_std_sort(Expr expr, OptimizedExpr)
{
    static_assert(std::is_same_v<decltype(stdSort<ProjExprFrequency, Compare>(expr)), OptimizedExpr>);
}

TEST(TestCTSmarts, StdSort)
{
    using C = AliphaticAtom<6>;
    using O = AliphaticAtom<8>;
    using S = AliphaticAtom<16>;

    auto less = ctll::list<S, O, C>{};
    auto greater = ctll::list<C, O, S>{};

    auto permutation1 = ctll::list<C, O, S>{};
    auto permutation2 = ctll::list<C, S, O>{};
    auto permutation3 = ctll::list<O, C, S>{};
    auto permutation4 = ctll::list<O, S, C>{};
    auto permutation5 = ctll::list<S, C, O>{};
    auto permutation6 = ctll::list<S, O, C>{};

    test_std_sort<std::less<>>(permutation1, less);
    test_std_sort<std::less<>>(permutation2, less);
    test_std_sort<std::less<>>(permutation3, less);
    test_std_sort<std::less<>>(permutation4, less);
    test_std_sort<std::less<>>(permutation5, less);
    test_std_sort<std::less<>>(permutation6, less);

    test_std_sort<std::greater<>>(permutation1, greater);
    test_std_sort<std::greater<>>(permutation2, greater);
    test_std_sort<std::greater<>>(permutation3, greater);
    test_std_sort<std::greater<>>(permutation4, greater);
    test_std_sort<std::greater<>>(permutation5, greater);
    test_std_sort<std::greater<>>(permutation6, greater);
}
*/

template<typename Goal = BestCase, typename Expr, typename OptimizedExpr>
constexpr void test_optimize_expression(Expr expr, OptimizedExpr)
{
    static_assert(std::is_same_v<decltype(optimizeExpression<Goal>(expr)), OptimizedExpr>);
}

TEST(TestCTSmarts, OptimizeExpression)
{
    // 6    0.49
    // 8    0.17
    // 16   0.00543654

    using C = AliphaticAtom<6>;
    using O = AliphaticAtom<8>;
    using S = AliphaticAtom<16>;

    // n = 1

    test_optimize_expression(C{}, C{});

    // n = 2

    auto less_2 = ctll::list<O, C>{};
    auto greater_2 = ctll::list<C, O>{};
    auto perm_2_1 = ctll::list<C, O>{};
    auto perm_2_2 = ctll::list<O, C>{};

    test_optimize_expression(And{perm_2_1}, And{less_2});
    test_optimize_expression(And{perm_2_1}, And{less_2});

    test_optimize_expression(Or{perm_2_1}, Or{greater_2});
    test_optimize_expression(Or{perm_2_2}, Or{greater_2});

    test_optimize_expression<WorstCase>(And{perm_2_1}, And{greater_2});
    test_optimize_expression<WorstCase>(And{perm_2_1}, And{greater_2});

    test_optimize_expression<WorstCase>(Or{perm_2_1}, Or{less_2});
    test_optimize_expression<WorstCase>(Or{perm_2_2}, Or{less_2});

    // n = 3

    auto less_3 = ctll::list<S, O, C>{};
    auto greater_3 = ctll::list<C, O, S>{};
    auto perm_3_1 = ctll::list<C, O, S>{};
    auto perm_3_2 = ctll::list<C, S, O>{};
    auto perm_3_3 = ctll::list<O, C, S>{};
    auto perm_3_4 = ctll::list<O, S, C>{};
    auto perm_3_5 = ctll::list<S, C, O>{};
    auto perm_3_6 = ctll::list<S, O, C>{};

    test_optimize_expression(And{perm_3_1}, And{less_3});
    test_optimize_expression(And{perm_3_2}, And{less_3});
    test_optimize_expression(And{perm_3_3}, And{less_3});
    test_optimize_expression(And{perm_3_4}, And{less_3});
    test_optimize_expression(And{perm_3_5}, And{less_3});
    test_optimize_expression(And{perm_3_6}, And{less_3});

    test_optimize_expression(Or{perm_3_1}, Or{greater_3});
    test_optimize_expression(Or{perm_3_2}, Or{greater_3});
    test_optimize_expression(Or{perm_3_3}, Or{greater_3});
    test_optimize_expression(Or{perm_3_4}, Or{greater_3});
    test_optimize_expression(Or{perm_3_5}, Or{greater_3});
    test_optimize_expression(Or{perm_3_6}, Or{greater_3});

    test_optimize_expression<WorstCase>(And{perm_3_1}, And{greater_3});
    test_optimize_expression<WorstCase>(And{perm_3_2}, And{greater_3});
    test_optimize_expression<WorstCase>(And{perm_3_3}, And{greater_3});
    test_optimize_expression<WorstCase>(And{perm_3_4}, And{greater_3});
    test_optimize_expression<WorstCase>(And{perm_3_5}, And{greater_3});
    test_optimize_expression<WorstCase>(And{perm_3_6}, And{greater_3});

    test_optimize_expression<WorstCase>(Or{perm_3_1}, Or{less_3});
    test_optimize_expression<WorstCase>(Or{perm_3_2}, Or{less_3});
    test_optimize_expression<WorstCase>(Or{perm_3_3}, Or{less_3});
    test_optimize_expression<WorstCase>(Or{perm_3_4}, Or{less_3});
    test_optimize_expression<WorstCase>(Or{perm_3_5}, Or{less_3});
    test_optimize_expression<WorstCase>(Or{perm_3_6}, Or{less_3});

    // Not

    auto best_not = Not<And<C, O>>{};
    auto worst_not = Not<And<O, C>>{};

    test_optimize_expression(best_not, best_not);
    test_optimize_expression(worst_not, best_not);

    test_optimize_expression<WorstCase>(best_not, worst_not);
    test_optimize_expression<WorstCase>(worst_not, worst_not);

    // depth > 1

    auto best_depth = And<Or<O, S>, C>{};
    auto worst_depth = And<C, Or<S, O>>{};

    static_assert(expressionFrequency(best_depth) == expressionFrequency(O{}));

    auto depth_perm_1 = And<C, Or<O, S>>{};
    auto depth_perm_2 = And<C, Or<S, O>>{};
    auto depth_perm_3 = And<Or<O, S>, C>{};
    auto depth_perm_4 = And<Or<S, O>, C>{};

    test_optimize_expression(depth_perm_1, best_depth);
    test_optimize_expression(depth_perm_2, best_depth);
    test_optimize_expression(depth_perm_3, best_depth);
    test_optimize_expression(depth_perm_4, best_depth);

    test_optimize_expression<WorstCase>(depth_perm_1, worst_depth);
    test_optimize_expression<WorstCase>(depth_perm_2, worst_depth);
    test_optimize_expression<WorstCase>(depth_perm_3, worst_depth);
    test_optimize_expression<WorstCase>(depth_perm_4, worst_depth);
}


template<typename Goal = BestCase, typename Expr, typename OptimizedExpr>
constexpr void test_optimize_expressions(Expr expr, OptimizedExpr)
{
    static_assert(std::is_same_v<decltype(optimizeExpressions<Goal>(expr)), OptimizedExpr>);
}

TEST(TestCTSmarts, OptimizeExpressions)
{
    // 6    0.49
    // 8    0.17
    // 16   0.00543654

    using C = AliphaticAtom<6>;
    using O = AliphaticAtom<8>;
    using S = AliphaticAtom<16>;

    using Best_1 = And<O, C>;
    using Worst_1 = And<C, O>;
    using Best_2 = And<S, C>;
    using Worst_2 = And<C, S>;
    using Best_3 = And<S, O>;
    using Worst_3 = And<O, S>;

    test_optimize_expression(Best_1{}, Best_1{});
    test_optimize_expression(Worst_1{}, Best_1{});
    test_optimize_expression(Best_2{}, Best_2{});
    test_optimize_expression(Worst_2{}, Best_2{});
    test_optimize_expression(Best_3{}, Best_3{});
    test_optimize_expression(Worst_3{}, Best_3{});

    auto perm1 = ctll::list<Best_1 , Best_2 , Best_3 >{};
    auto perm2 = ctll::list<Best_1 , Best_2 , Worst_3>{};
    auto perm3 = ctll::list<Best_1 , Worst_2, Best_3 >{};
    auto perm4 = ctll::list<Best_1 , Worst_2, Worst_3>{};
    auto perm5 = ctll::list<Worst_1, Best_2 , Best_3 >{};
    auto perm6 = ctll::list<Worst_1, Best_2 , Worst_3>{};
    auto perm7 = ctll::list<Worst_1, Worst_2, Best_3 >{};
    auto perm8 = ctll::list<Worst_1, Worst_2 , Worst_3>{};

    auto best = perm1;
    auto worst = perm8;

    test_optimize_expressions(perm1, best);
    test_optimize_expressions(perm2, best);
    test_optimize_expressions(perm3, best);
    test_optimize_expressions(perm4, best);
    test_optimize_expressions(perm5, best);
    test_optimize_expressions(perm6, best);
    test_optimize_expressions(perm7, best);
    test_optimize_expressions(perm8, best);

    test_optimize_expressions<WorstCase>(perm1, worst);
    test_optimize_expressions<WorstCase>(perm2, worst);
    test_optimize_expressions<WorstCase>(perm3, worst);
    test_optimize_expressions<WorstCase>(perm4, worst);
    test_optimize_expressions<WorstCase>(perm5, worst);
    test_optimize_expressions<WorstCase>(perm6, worst);
    test_optimize_expressions<WorstCase>(perm7, worst);
    test_optimize_expressions<WorstCase>(perm8, worst);
}


TEST(TestCTSmarts, AtomFrequency)
{
    using C = AliphaticAtom<6>;
    using N = AliphaticAtom<7>;
    using O = AliphaticAtom<8>;

    auto smarts = Smarts<"CNO">{};

    constexpr auto atomFreq = AtomFrequency{smarts}.data;

    static_assert(atomFreq.size() == 3);
    static_assert(atomFreq[0] == expressionFrequency(C{}));
    static_assert(atomFreq[1] == expressionFrequency(N{}));
    static_assert(atomFreq[2] == expressionFrequency(O{}));
}

TEST(TestCTSmarts, IncidentList)
{
    /*
    using C = AliphaticAtom<6>;
    using N = AliphaticAtom<7>;
    using O = AliphaticAtom<8>;
    using F = AliphaticAtom<9>;
    */

    auto smarts = Smarts<"CN(O)F">{};

    constexpr auto edgeList = EdgeList{smarts};
    constexpr auto vertexDegree = VertexDegree{smarts, edgeList};
    constexpr auto incidentList = IncidentList{smarts, edgeList, vertexDegree};

    static_assert(incidentList.data.size() == 4 * 3);
    static_assert(incidentList.data[0] == 0);
    static_assert(incidentList.data[3] == 0);
    static_assert(incidentList.data[4] == 1);
    static_assert(incidentList.data[5] == 2);
    static_assert(incidentList.data[6] == 1);
    static_assert(incidentList.data[9] == 2);


    constexpr auto atomFreq = AtomFrequency{smarts};
    constexpr auto optimizeIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq};

    static_assert(optimizeIncidentList.data.size() == 4 * 3);
    static_assert(optimizeIncidentList.data[0] == 0);
    static_assert(optimizeIncidentList.data[3] == 2);
    static_assert(optimizeIncidentList.data[4] == 1);
    static_assert(optimizeIncidentList.data[5] == 0);
    static_assert(optimizeIncidentList.data[6] == 1);
    static_assert(optimizeIncidentList.data[9] == 2);

}


TEST(TestCTSmarts, Real)
{
    static_assert(Real<1>::value == 1.0);
    static_assert(Real<1, 0>::value == 1.0);
    static_assert(Real<1, 1>::value == 10.0);
    static_assert(Real<1, -1>::value == 0.1);
    static_assert(Real<2, 2>::value == 200.0);
    static_assert(Real<2, -2>::value == 0.02);
}

TEST(TestCTSmarts, Merge)
{
    using A = Number<0>;
    using B = Number<1>;
    using C = Number<2>;
    using D = Number<3>;

    static_assert(std::is_same_v<decltype(merge(ctll::list<A, B>{}, ctll::list<C, D>{})), ctll::list<A, B, C, D>>);
    static_assert(std::is_same_v<decltype(merge(ctll::list<A, B>{}, ctll::list<A, B>{})), ctll::list<A, B>>);
    static_assert(std::is_same_v<decltype(merge(ctll::list<A, C, B>{}, ctll::list<B, D, C>{})), ctll::list<A, B, D, C>>);


    //identify_type<decltype(merge(ctll::list<A, C, B>{}, ctll::list<B, D, C>{}))>{};

    //static_assert(std::is_same_v<decltype(merge(ctll::list<bool, int>{}, ctll::list<float, double>{})), ctll::list<bool, int, float, double>>);


}

template<ctll::fixed_string SMARTS, typename Primitives>
constexpr auto test_required_atom_primitives() noexcept
{
    auto smarts = Smarts<SMARTS>{};
    static_assert(std::is_same_v<decltype(requiredAtomPrimitives(smarts.atoms)), Primitives>);

}

TEST(TestCTSmarts, RequiredAtomPrimitives)
{
    using C = AliphaticAtom<6>;
    using N = AliphaticAtom<7>;
    using O = AliphaticAtom<8>;
    using F = AliphaticAtom<9>;

    auto smarts1 = Smarts<"CCO">{};

    test_required_atom_primitives<"C", ctll::list<C>>();
    test_required_atom_primitives<"N", ctll::list<N>>();
    test_required_atom_primitives<"O", ctll::list<O>>();
    test_required_atom_primitives<"F", ctll::list<F>>();
    test_required_atom_primitives<"CCCC", ctll::list<C>>();
    test_required_atom_primitives<"CCO", ctll::list<C, O>>();
    test_required_atom_primitives<"CCN", ctll::list<C, N>>();
    test_required_atom_primitives<"CNOF", ctll::list<C, N, O, F>>();
    test_required_atom_primitives<"CNOFCNOF", ctll::list<C, N, O, F>>();

}



TEST(TestCTSmarts, NumAtomBondFilter)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]
    auto smarts1 = Smarts<"CCC">{};
    auto smarts2 = Smarts<"CCCCC">{}; // more atoms than mol

    auto filterPolicy = FilterPolicy<NumAtomBondFilter>{};
    auto filterPolicyHelper1 = FilterPolicyHelper{smarts1, filterPolicy.filters};
    auto filterPolicyHelper2 = FilterPolicyHelper{smarts2, filterPolicy.filters};

    EXPECT_FALSE(filterPolicyHelper1.reject(mol));
    EXPECT_TRUE(filterPolicyHelper2.reject(mol));
}

TEST(TestCTSmarts, ElementFilter)
{

    auto mol = mockAcetateAnion(); // CC(=O)[O-]
    auto smarts1 = Smarts<"CCO">{};
    auto smarts2 = Smarts<"CCN">{};
    auto smarts3 = Smarts<"[Ti]">{};

    auto allElements1 = impl::enabledElements<Real<1>>(smarts1);
    auto allElements2 = impl::enabledElements<Real<1>>(smarts2);
    auto allElements3 = impl::enabledElements<Real<1>>(smarts3);

    static_assert(std::is_same_v<decltype(allElements1), ctll::list<Number<6>, Number<8>>>);
    static_assert(std::is_same_v<decltype(allElements2), ctll::list<Number<6>, Number<7>>>);
    static_assert(std::is_same_v<decltype(allElements3), ctll::list<Number<22>>>);

    auto rareElements1 = impl::enabledElements<Real<1, -2>>(smarts1);
    auto rareElements2 = impl::enabledElements<Real<1, -2>>(smarts2);
    auto rareElements3 = impl::enabledElements<Real<1, -2>>(smarts3);

    static_assert(std::is_same_v<decltype(rareElements1), ctll::empty_list>);
    static_assert(std::is_same_v<decltype(rareElements2), ctll::empty_list>);
    static_assert(std::is_same_v<decltype(rareElements3), ctll::list<Number<22>>>);


    auto filterPolicy = FilterPolicy<ElementFilter<Real<1>>>{};
    auto filterPolicyHelper1 = FilterPolicyHelper{smarts1, filterPolicy.filters};
    auto filterPolicyHelper2 = FilterPolicyHelper{smarts2, filterPolicy.filters};
    auto filterPolicyHelper3 = FilterPolicyHelper{smarts3, filterPolicy.filters};

    EXPECT_FALSE(filterPolicyHelper1.reject(mol));
    EXPECT_TRUE(filterPolicyHelper2.reject(mol));
    EXPECT_TRUE(filterPolicyHelper3.reject(mol));

}





TEST(TestCTSmarts, Debug)
{
    auto mol = mockAcetateAnion(); // CC(=O)[O-]

    auto smarts1 = Smarts<"CCC">{};
    auto smarts2 = Smarts<"CCCCC">{};

    auto filterPolicy = FilterPolicy<NumAtomBondFilter>{};

    auto filterPolicyHelper1 = FilterPolicyHelper{smarts1, filterPolicy.filters};
    auto filterPolicyHelper2 = FilterPolicyHelper{smarts2, filterPolicy.filters};

    EXPECT_FALSE(filterPolicyHelper1.reject(mol));
    EXPECT_TRUE(filterPolicyHelper2.reject(mol));

    auto filters1 = filterPolicyHelper1.filters;

}













void test_Errors()
{
    test_parse_error<"CC[]", EmptyBracketAtomError<Pos<3>> >(); // FIXME: pos 2 = start
    //SmartsT<"[]"> a;


}





#include <iostream>


#include <tuple>
#include <vector>
#include <array>









//#include <ctsmarts/FastMolecule.hpp>

//using namespace std;









template <ctll::fixed_string input>
constexpr bool test()
{
    return ctll::parser<SmartsGrammar, input, SmartsActions>::template correct_with<>;
}






















template <ctll::fixed_string input, typename Atom, bool Expected = true>
constexpr void test_atom_match()
{
    using AST = ctll::parser<SmartsGrammar, input, SmartsActions>::template output<SmartsContext<>>;
    static_assert(AST::is_correct);
    using smarts = AST::output_type;
    static_assert(ctll::size(smarts::atoms) == 1);
    auto expr = ctll::front(smarts::atoms);
    constexpr auto m = contains(Atom(), expr);
    static_assert(m == Expected);
}































//template<typename T> /*requires
//{ ctll::size(std::list<>())+0 } -> std::same_as<std::size_t>; }

template<typename Molecule, typename Atom, typename ...Ts>
auto matchAtomExpr(const Molecule &mol, const Atom &atom, ctll::list<Ts...> atoms, std::size_t i)
{
    return with_n<ctll::size(atoms), bool>(i, [&mol, &atom, atoms] (auto j) {
        return matchAtomExpr(mol, atom, get<j>(atoms));
    });
}

template<typename T>
struct identify_type;


/*
template<ctll::fixed_string Str, auto I, auto N>
constexpr auto toCharArray2(Str, std::array<char, N> a)
{
    if constexpr (I == N)
        return a;
    else {
        //array[I] =
        return toCharArray2<Str, I + 1>(str, a);
    }
}
*/

template<ctll::fixed_string Str>
constexpr auto toCharArray()
{
    auto a = std::array<char, Str.size()>{};
    //return toCharArray2<Str, 0>(str, a);
    return a;
}


/*
template<typename Str>
constexpr auto toCharArray(Str)
{
    return std::array<char, Str().size()>



}
*/

/*
template<ctll::fixed_string Str>
struct Msg
{
    static constexpr inline std::

};
*/



template<typename Molecule>
constexpr auto carboxylCarbon = [] (Molecule &mol) { return captures<"C(=O)O", Molecule>(mol); };





template<typename T, typename U>
std::ostream& operator<<(std::ostream &os, const std::pair<T, U> &p)
{
    os << "( " << p.first << ", " << p.second << " )";
    return os;
}

/* Test
std::ostream& operator<<(std::ostream &os, const FM::Offset &offset)
{
    os << "Offset( valid = " << offset.valid
       << ", value = " << offset.value
       << ", n = " << offset.n
       << ", m = " << offset.m
       << ", i = " << offset.i
       << " )";
    return os;
}
*/


    /* TEST
    std::cout << FM::offset(FM::AtomArray(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::AtomStruct(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::BondArray(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::BondStruct(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::IncidentStruct(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::IncidentIndex(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::IncidentArray(), FM::MoleculeStruct()) << std::endl;

    std::cout << FM::sizeOf(FM::MoleculeStruct()) << std::endl;
    std::cout << FM::sizeOf(FM::MoleculeStruct(), 0, 0) << std::endl;
    std::cout << FM::sizeOf(FM::MoleculeStruct(), 1, 0) << std::endl;
    std::cout << FM::sizeOf(FM::MoleculeStruct(), 0, 1) << std::endl;
    std::cout << FM::sizeOf(FM::MoleculeStruct(), 1, 1) << std::endl;
    std::cout << FM::sizeOf(FM::MoleculeStruct(), 2, 2) << std::endl;

    std::cout << FM::offset(FM::ElementData(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::SourceData(), FM::MoleculeStruct()) << std::endl;
    std::cout << FM::offset(FM::TargetData(), FM::MoleculeStruct()) << std::endl;

    std::cout << FM::offset(FM::ElementData(), FM::MoleculeStruct(), 1, 0, 0) << std::endl;
    std::cout << FM::offset(FM::IsotopeData(), FM::MoleculeStruct(), 1, 0, 0) << std::endl;

    std::cout << FM::offset(FM::ElementData(), FM::MoleculeStruct(), 2, 0, 1) << std::endl;
    std::cout << FM::offset(FM::IsotopeData(), FM::MoleculeStruct(), 2, 0, 1) << std::endl;
    */



















    /* TEST
    auto [match, N, O, C] = capture<"[N:1]([C:3])[O:2]">(mol);

    std::cout << C->GetType() << std::endl;
    std::cout << N->GetType() << std::endl;
    std::cout << O->GetType() << std::endl;
*/


    //identify_type<decltype(atoms)>();

    /*
    constexpr auto cap = captureMapping(smarts);

    std::cout << "capute map:" << std::endl;
    for (auto i = 0; i < cap.size(); ++i)
        std::cout << "    " << i << " -> " << cap[i] << std::endl;

    Isomorphism iso(mol, smarts);

    auto map = iso.single();

    std::cout << "isomorphism map:" << std::endl;
        for (auto i = 0; i < map.size(); ++i)
        std::cout << "    " << i << " -> " << map[i] << std::endl;

    auto N = get_atom(mol, map[cap[0]]);
    auto O = get_atom(mol, map[cap[1]]);
    auto C = get_atom(mol, map[cap[2]]);

    std::cout << C->GetType() << std::endl;
    std::cout << N->GetType() << std::endl;
    std::cout << O->GetType() << std::endl;
*/






    /*

    auto smarts = Smarts<"[N:1]([C:1])[O:3]">();

    Isomorphism iso(mol, smarts);

    std::cout << iso.contains() << std::endl;;
*/



    /*
    auto smarts = Smarts<"CC">();


    auto input = smarts.input();
    auto context = smarts.context;
    auto valid = smarts.valid;
    auto error = smarts.error;
    auto position = smarts.position;




    MockMolecule mol;
    auto C0 = mol.addAtom();
    auto C1 = mol.addAtom();
    //auto C2 = mol.addAtom();
    mol.addBond(0, 1);
    //mol.addBond(1, 2);
    //mol.addBond(0, 2);





    NoMapping mapping1;
    Isomorphism iso1(mol, smarts);
    iso1.match(mapping1);

    std::cout << mapping1.match << std::endl;

    NoMapping mapping2;
    IsomorphismRecursive iso2(mol, smarts);
    iso2.match(mapping2);

    std::cout << mapping2.match << std::endl;

*/



/*
    test_atom_match<"C", MockAtom<6>>();
    test_atom_match<"C", MockAtom<7>, false>();


    test_atom_match<"[N,O]", MockAtom<7>>();
    test_atom_match<"[N,O]", MockAtom<8>>();
    test_atom_match<"[N,O]", MockAtom<6>, false>();

    test_atom_match<"[N,O,S]", MockAtom<7>>();
    test_atom_match<"[N,O,S]", MockAtom<8>>();
    test_atom_match<"[N,O,S]", MockAtom<16>>();
    test_atom_match<"[N,O,S]", MockAtom<6>, false>();
*/





    //static_assert(std::same_as<smarts, void>);



    //static_assert(test<"C">());




