#include "TestCTSmarts.hpp"

#ifdef KITIMAR_WITH_OPENBABEL
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <openbabel/parsmart.h>
#endif



using namespace Kitimar;
using namespace Kitimar::CTSmarts;



template<typename T>
struct identify_type;

TEST_CASE("EdgeList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    constexpr auto edges = EdgeList{smarts}.data;

    static_assert(edges.size() == 4);
    static_assert(edges[0] == Edge{0, 0, 1});
    static_assert(edges[1] == Edge{1, 1, 2});
    static_assert(edges[2] == Edge{2, 1, 3});
    static_assert(edges[3] == Edge{3, 3, 0});
}


#ifdef KITIMAR_WITH_IOSTREAM

TEST_CASE("DfsSearch")
{
    auto smarts = Smarts<"*1**12**2*">{};

    auto edges = EdgeList(smarts);
    auto degrees = VertexDegree(smarts, edges);
    auto adjList = IncidentList(smarts, edges, degrees);

    auto dfsSearchEvents = DfsSearchEvents(smarts, adjList);

    for (const auto &event : dfsSearchEvents.events)
        std::cout << event << std::endl;
}

#endif // KITIMAR_WITH_IOSTREAM


TEST_CASE("ExpressionFrequency")
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

TEST_CASE("CtllSort")
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

TEST_CASE("StdSort)
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

TEST_CASE("OptimizeExpression")
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

TEST_CASE("OptimizeExpressions")
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


TEST_CASE("AtomFrequency")
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

TEST_CASE("IncidentList")
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


TEST_CASE("Real")
{
    static_assert(Real<1>::value == 1.0);
    static_assert(Real<1, 0>::value == 1.0);
    static_assert(Real<1, 1>::value == 10.0);
    static_assert(Real<1, -1>::value == 0.1);
    static_assert(Real<2, 2>::value == 200.0);
    static_assert(Real<2, -2>::value == 0.02);
}

TEST_CASE("Merge")
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

TEST_CASE("RequiredAtomPrimitives")
{
    using C = AliphaticAtom<6>;
    using N = AliphaticAtom<7>;
    using O = AliphaticAtom<8>;
    using F = AliphaticAtom<9>;

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



TEST_CASE("NumAtomBondFilter")
{
    auto mol = Molecule::mockAcetateAnion(); // CC(=O)[O-]
    auto smarts1 = Smarts<"CCC">{};
    auto smarts2 = Smarts<"CCCCC">{}; // more atoms than mol

    auto filterPolicy = FilterPolicy<NumAtomBondFilter>{};
    auto filterPolicyHelper1 = FilterPolicyHelper{smarts1, filterPolicy.filters};
    auto filterPolicyHelper2 = FilterPolicyHelper{smarts2, filterPolicy.filters};

    CHECK(!filterPolicyHelper1.reject(mol));
    CHECK(filterPolicyHelper2.reject(mol));
}

TEST_CASE("ElementFilter")
{
    auto mol = Molecule::mockAcetateAnion(); // CC(=O)[O-]
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

    CHECK(!filterPolicyHelper1.reject(mol));
    CHECK(filterPolicyHelper2.reject(mol));
    CHECK(filterPolicyHelper3.reject(mol));

}


#include <catch2/benchmark/catch_benchmark.hpp>



/*
constexpr auto benchmarkSMARTS = std::make_tuple(
    ctll::fixed_string{"CC"},
    ctll::fixed_string{"CO"}
);

using BenchmarkSMARTS = std::remove_cvref_t<decltype(benchmarkSMARTS)>;
*/

template<ctll::fixed_string SMARTS>
struct FixedString
{
        static constexpr inline auto smarts = SMARTS;
};

using BenchmarkSMARTS = std::tuple<
    FixedString<"CC">,
    FixedString<"CO">
>;


TEMPLATE_LIST_TEST_CASE("Benchmark", "", BenchmarkSMARTS)
{
    auto SMARTS = TestType::smarts;

    BENCHMARK(Util::toString(SMARTS)) {
        return 1;
        //std::cout << Util::toString(SMARTS) << std::endl;
    };

}

TEST_CASE("Debug")
{


    /*
    for (auto i = 0; i < 10; ++i) {
        BENCHMARK("foo") {
            return 1;
        };
    }
    */







}











