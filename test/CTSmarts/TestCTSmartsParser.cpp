#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

template <ctll::fixed_string SMARTS, typename AtomExpr/*, int AtomIndex = 0*/>
constexpr void test_atom_expr()
{
    constexpr auto AtomIndex = 0;
    auto smarts = Smarts<SMARTS>{};
    static_assert(ctll::size(smarts.atoms) > AtomIndex);
    auto expr = get<AtomIndex>(smarts.atoms);
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

TEST_CASE("AliphaticAtom")
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

TEST_CASE("AromaticAtom")
{
    // 'b' | 'c' | 'n' | 'o' | 's' | 'p'
    test_atom_expr<"b",  AromaticAtom<5>>();
    test_atom_expr<"c",  AromaticAtom<6>>();
    test_atom_expr<"n",  AromaticAtom<7>>();
    test_atom_expr<"o",  AromaticAtom<8>>();
    test_atom_expr<"p",  AromaticAtom<15>>();
    test_atom_expr<"s",  AromaticAtom<16>>();
}

TEST_CASE("AnyAtom")
{
    test_atom_expr<"*",  AnyAtom>();
    test_atom_expr<"a",  AnyAromatic>();
    test_atom_expr<"A",  AnyAliphatic>();
}

TEST_CASE("BracketAliphaticAtom")
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

TEST_CASE("BracketAromaticAtom")
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

TEST_CASE("BracketAnyAtom")
{
    test_atom_expr<"[*]", AnyAtom >();
    test_atom_expr<"[a]", AnyAromatic >();
    test_atom_expr<"[A]", AnyAliphatic >();
}

TEST_CASE("Isotope")
{
    test_atom_expr<"[1]",   Isotope<1> >();
    test_atom_expr<"[2]",   Isotope<2> >();
    test_atom_expr<"[0]",   Isotope<0> >();
    test_atom_expr<"[42]",  Isotope<42> >();
    test_atom_expr<"[123]", Isotope<123> >();
}

TEST_CASE("Element")
{
    test_atom_expr<"[#1]",   Element<1> >();
    test_atom_expr<"[#2]",   Element<2> >();
    test_atom_expr<"[#0]",   Element<0> >();
    test_atom_expr<"[#42]",  Element<42> >();
    test_atom_expr<"[#123]", Element<123> >();
}

TEST_CASE("Degree")
{
    test_atom_expr<"[D]",   Degree<1> >();
    test_atom_expr<"[D1]",  Degree<1> >();
    test_atom_expr<"[D2]",  Degree<2> >();
    test_atom_expr<"[D0]",  Degree<0> >();
    test_atom_expr<"[D42]", Degree<42> >();
}

TEST_CASE("Valence")
{
    test_atom_expr<"[v]",   Valence<1> >();
    test_atom_expr<"[v1]",  Valence<1> >();
    test_atom_expr<"[v2]",  Valence<2> >();
    test_atom_expr<"[v0]",  Valence<0> >();
    test_atom_expr<"[v42]", Valence<42> >();
}

TEST_CASE("Connectivity")
{
    test_atom_expr<"[X]",   Connectivity<1> >();
    test_atom_expr<"[X1]",  Connectivity<1> >();
    test_atom_expr<"[X2]",  Connectivity<2> >();
    test_atom_expr<"[X0]",  Connectivity<0> >();
    test_atom_expr<"[X42]", Connectivity<42> >();
}

TEST_CASE("TotalH")
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

TEST_CASE("ImplicitH")
{
    test_atom_expr<"[h]",  ImplicitH<1> >();
    test_atom_expr<"[h0]", ImplicitH<0> >();
    test_atom_expr<"[h1]", ImplicitH<1> >();
    test_atom_expr<"[h3]", ImplicitH<3> >();
}

TEST_CASE("CyclicAtom")
{
    test_atom_expr<"[R]", Cyclic >();
    test_atom_expr<"[r]", Cyclic >();
    test_atom_expr<"[x]", Cyclic >();
}

TEST_CASE("AcyclicAtom")
{
    //test_atom_expr<"[R0]", Acyclic >();  FIXME
    test_atom_expr<"[r0]", Acyclic >();
    test_atom_expr<"[x0]", Acyclic >();
}

TEST_CASE("RingCount")
{
    test_atom_expr<"[R1]",  RingCount<1> >();
    test_atom_expr<"[R2]",  RingCount<2> >();
    test_atom_expr<"[R42]", RingCount<42> >();
}

TEST_CASE("RingSize")
{
    test_atom_expr<"[r1]",  RingSize<1> >();
    test_atom_expr<"[r2]",  RingSize<2> >();
    test_atom_expr<"[r42]", RingSize<42> >();
}

TEST_CASE("RingConnectivity")
{
    test_atom_expr<"[x1]",  RingConnectivity<1> >();
    test_atom_expr<"[x2]",  RingConnectivity<2> >();
    test_atom_expr<"[x42]", RingConnectivity<42> >();
}


TEST_CASE("Charge")
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

TEST_CASE("Chiral")
{
}

TEST_CASE("Class")
{
}

TEST_CASE("BracketAtom")
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


TEST_CASE("BondPrimitive")
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

TEST_CASE("BondExpr")
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

TEST_CASE("RingBond")
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

TEST_CASE("Operators")
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

template<typename T>
struct identify_type;

template<typename Expr>
using R1 = BasicSmarts<ctll::list<Atom<0, Expr>>>;

template<typename Expr1, typename Expr2>
using R2 = BasicSmarts<ctll::list<Atom<0, Expr1>, Atom<1, Expr2>>, ctll::list<Bond<0, 0, 1, ImplicitBond>>>;

TEST_CASE("Recursive")
{
    using S = AliphaticAtom<16>;
    using N = AliphaticAtom<7>;

    test_atom_expr<"[S]", S >();
    test_atom_expr<"[S$(*)]", And<S, R1<AnyAtom>> >();
    test_atom_expr<"[S$(*N)]", And<S, R2<AnyAtom, N>> >();
    test_atom_expr<"[S!$(*N)]", And<S, Not<R2<AnyAtom, N>>> >();


}


TEST_CASE("GetCentralAtom")
{
    CHECK(getCentralAtom(Smarts<"C">{}) == -1);
    CHECK(getCentralAtom(Smarts<"CC">{}) == -1);
    CHECK(getCentralAtom(Smarts<"CCC">{}) == 1);
    CHECK(getCentralAtom(Smarts<"C(C)C">{}) == 0);
    CHECK(getCentralAtom(Smarts<"CCCC">{}) == -1);
    CHECK(getCentralAtom(Smarts<"CC(C)C">{}) == 1);
    CHECK(getCentralAtom(Smarts<"C(C)(C)C">{}) == 0);
    CHECK(getCentralAtom(Smarts<"CCCCC">{}) == -1);
    CHECK(getCentralAtom(Smarts<"CC(C)CC">{}) == -1);
    CHECK(getCentralAtom(Smarts<"CC(C)(C)C">{}) == 1);
    CHECK(getCentralAtom(Smarts<"C(C)(C)(C)C">{}) == 0);
    CHECK(getCentralAtom(Smarts<"C1CC1">{}) == -1);
}
