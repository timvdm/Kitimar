#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../Test.hpp"

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

//
//
//

TEST_CASE("AtomFequency")
{
    auto smarts = Smarts<"CN(O)F">{};

    constexpr auto freq = AtomFrequency{smarts}.data;

    static_assert(freq.size() == 4);
    static_assert(freq[0] == expressionFrequency(AliphaticAtom<6>{}));
    static_assert(freq[1] == expressionFrequency(AliphaticAtom<7>{}));
    static_assert(freq[2] == expressionFrequency(AliphaticAtom<8>{}));
    static_assert(freq[3] == expressionFrequency(AliphaticAtom<9>{}));
}

//
// OptimizeIncidentList
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueOptimizeIncidentList")
{
    auto smarts = Smarts<"CN(O)F">{};

    constexpr auto edgeList = ValueEdgeList{smarts};
    constexpr auto vertexDegree = ValueVertexDegree{smarts, edgeList};
    constexpr auto incidentList = ValueIncidentList{smarts, edgeList, vertexDegree};
    constexpr auto atomFreq = AtomFrequency{smarts};
    constexpr auto optimizeIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq, std::false_type{}};

    static_assert(optimizeIncidentList.data.size() == 4 * 3);
    static_assert(optimizeIncidentList.data[0] == 0);
    static_assert(optimizeIncidentList.data[3] == 2);
    static_assert(optimizeIncidentList.data[4] == 1);
    static_assert(optimizeIncidentList.data[5] == 0);
    static_assert(optimizeIncidentList.data[6] == 1);
    static_assert(optimizeIncidentList.data[9] == 2);
}

TEST_CASE("OptimizeIncidentList<SeedBond>")
{
    auto smarts = Smarts<"C(C)S">{};

    constexpr auto edgeList = ValueEdgeList{smarts};
    constexpr auto vertexDegree = ValueVertexDegree{smarts, edgeList};
    constexpr auto incidentList = ValueIncidentList{smarts, edgeList, vertexDegree};
    constexpr auto atomFreq = AtomFrequency{smarts};

    // no seed bond
    constexpr auto optimizeIncidentList1 = OptimizeIncidentList{smarts, incidentList, atomFreq, std::false_type{}};
    static_assert(optimizeIncidentList1.data.size() == 3 * 2);
    static_assert(optimizeIncidentList1.data[0] == 1);
    static_assert(optimizeIncidentList1.data[1] == 0);
    static_assert(optimizeIncidentList1.data[2] == 0);
    static_assert(optimizeIncidentList1.data[4] == 1);

    // seed bond
    constexpr auto optimizeIncidentList2 = OptimizeIncidentList{smarts, incidentList, atomFreq, std::true_type{}};
    static_assert(optimizeIncidentList2.data.size() == 3 * 2);
    static_assert(optimizeIncidentList2.data[0] == 0);
    static_assert(optimizeIncidentList2.data[1] == 1);
    static_assert(optimizeIncidentList2.data[2] == 0);
    static_assert(optimizeIncidentList2.data[4] == 1);
}

#endif // KITIMAR_VALUE_BASED
