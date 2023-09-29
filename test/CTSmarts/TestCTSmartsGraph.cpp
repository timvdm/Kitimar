#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../Test.hpp"

using namespace Kitimar;
using namespace Kitimar::CTSmarts;

//
// EdgeList
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueEdgeList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    constexpr auto edges = ValueEdgeList{smarts}.data;

    static_assert(edges.size() == 4);
    static_assert(edges[0] == ValueEdge{0, 0, 1});
    static_assert(edges[1] == ValueEdge{1, 1, 2});
    static_assert(edges[2] == ValueEdge{2, 1, 3});
    static_assert(edges[3] == ValueEdge{3, 3, 0});
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeEdgeList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    constexpr auto edges = TypeEdgeList{smarts}.data;

    static_assert(ctll::size(edges) == 4);
    static_assert(get<0>(edges).index == 0);
    static_assert(get<0>(edges).source == 0);
    static_assert(get<0>(edges).target == 1);
    static_assert(get<1>(edges).index == 1);
    static_assert(get<1>(edges).source == 1);
    static_assert(get<1>(edges).target == 2);
    static_assert(get<2>(edges).index == 2);
    static_assert(get<2>(edges).source == 1);
    static_assert(get<2>(edges).target == 3);
    static_assert(get<3>(edges).index == 3);
    static_assert(get<3>(edges).source == 3);
    static_assert(get<3>(edges).target == 0);
}

//
// VertexDegree
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueVertexDegree")
{
    auto smarts = Smarts<"**(*)**">{};

    constexpr auto edgeList = ValueEdgeList{smarts};
    constexpr auto degrees = ValueVertexDegree{smarts, edgeList}.data;

    static_assert(degrees.size() == 5);
    static_assert(degrees[0] == 1);
    static_assert(degrees[1] == 3);
    static_assert(degrees[2] == 1);
    static_assert(degrees[3] == 2);
    static_assert(degrees[4] == 1);
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeVertexDegree")
{
    auto smarts = Smarts<"**(*)**">{};

    constexpr auto edgeList = TypeEdgeList{smarts};
    constexpr auto degrees = TypeVertexDegree{smarts, edgeList}.data;

    static_assert(degrees.size() == 5);
    static_assert(degrees[0] == 1);
    static_assert(degrees[1] == 3);
    static_assert(degrees[2] == 1);
    static_assert(degrees[3] == 2);
    static_assert(degrees[4] == 1);
}

//
// IncidentList
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueIncidentList")
{
    auto smarts = Smarts<"CN(O)F">{};

    constexpr auto edgeList = ValueEdgeList{smarts};
    constexpr auto vertexDegree = ValueVertexDegree{smarts, edgeList};
    constexpr auto incidentList = ValueIncidentList{smarts, edgeList, vertexDegree};

    static_assert(incidentList.data.size() == 4 * 3);
    static_assert(incidentList.data[ 0] == 0);
    static_assert(incidentList.data[ 1] == -1);
    static_assert(incidentList.data[ 2] == -1);
    static_assert(incidentList.data[ 3] == 0);
    static_assert(incidentList.data[ 4] == 1);
    static_assert(incidentList.data[ 5] == 2);
    static_assert(incidentList.data[ 6] == 1);
    static_assert(incidentList.data[ 7] == -1);
    static_assert(incidentList.data[ 8] == -1);
    static_assert(incidentList.data[ 9] == 2);
    static_assert(incidentList.data[10] == -1);
    static_assert(incidentList.data[11] == -1);
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeIncidentList")
{
    auto smarts = Smarts<"CN(O)F">{};

    constexpr auto edgeList = TypeEdgeList{smarts};
    constexpr auto vertexDegree = TypeVertexDegree{smarts, edgeList};
    constexpr auto incidentList = TypeIncidentList{smarts, edgeList, vertexDegree};

    static_assert(incidentList.data.size() == 4 * 3);
    static_assert(incidentList.data[ 0] == 0);
    static_assert(incidentList.data[ 1] == -1);
    static_assert(incidentList.data[ 2] == -1);
    static_assert(incidentList.data[ 3] == 0);
    static_assert(incidentList.data[ 4] == 1);
    static_assert(incidentList.data[ 5] == 2);
    static_assert(incidentList.data[ 6] == 1);
    static_assert(incidentList.data[ 7] == -1);
    static_assert(incidentList.data[ 8] == -1);
    static_assert(incidentList.data[ 9] == 2);
    static_assert(incidentList.data[10] == -1);
    static_assert(incidentList.data[11] == -1);
}

//
// DfsSearch
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueDfsSearch")
{
    auto smarts = Smarts<"*1**12**2*">{};

    auto edges = ValueEdgeList{smarts};
    auto degrees = ValueVertexDegree{smarts, edges};
    auto incidentList = ValueIncidentList{smarts, edges, degrees};

    constexpr auto e = ValueDfsSearchEvents{smarts, incidentList}.events;

    static_assert(e[ 0] == ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, 0, true});
    static_assert(e[ 1] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 0, false});
    static_assert(e[ 2] == ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, 1, false});
    static_assert(e[ 3] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 1, false});
    static_assert(e[ 4] == ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, 2, false});
    static_assert(e[ 5] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 2, true});
    static_assert(e[ 6] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 2, true});
    static_assert(e[ 7] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 3, false});
    static_assert(e[ 8] == ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, 3, false});
    static_assert(e[ 9] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 4, false});
    static_assert(e[10] == ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, 4, false});
    static_assert(e[11] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 5, true});
    static_assert(e[12] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 5, true});
    static_assert(e[13] == ValueDfsSearchEvent{DfsSearchEventType::VisitEdge, 6, false});
    static_assert(e[14] == ValueDfsSearchEvent{DfsSearchEventType::VisitVertex, 5, false});
    static_assert(e[15] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, 5, false});
    static_assert(e[16] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 6, false});
    static_assert(e[17] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, 4, false});
    static_assert(e[18] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 4, false});
    static_assert(e[19] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, 3, false});
    static_assert(e[20] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 3, false});
    static_assert(e[21] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, 2, false});
    static_assert(e[22] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 1, false});
    static_assert(e[23] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, 1, false});
    static_assert(e[24] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackEdge, 0, false});
    static_assert(e[25] == ValueDfsSearchEvent{DfsSearchEventType::BacktrackVertex, 0, true});
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeDfsSearch")
{
    auto smarts = Smarts<"*1**12**2*">{};

    auto edges = TypeEdgeList{smarts};
    auto degrees = TypeVertexDegree{smarts, edges};
    auto incidentList = TypeIncidentList{smarts, edges, degrees};

    constexpr auto e = TypeDfsSearchEvents{smarts, incidentList}.events;

    static_assert(std::is_same_v< decltype(get< 0>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, 0, true> >);
    static_assert(std::is_same_v< decltype(get< 1>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 0, false> >);
    static_assert(std::is_same_v< decltype(get< 2>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, 1, false> >);
    static_assert(std::is_same_v< decltype(get< 3>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 1, false> >);
    static_assert(std::is_same_v< decltype(get< 4>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, 2, false> >);
    static_assert(std::is_same_v< decltype(get< 5>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 2, true> >);
    static_assert(std::is_same_v< decltype(get< 6>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 2, true> >);
    static_assert(std::is_same_v< decltype(get< 7>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 3, false> >);
    static_assert(std::is_same_v< decltype(get< 8>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, 3, false> >);
    static_assert(std::is_same_v< decltype(get< 9>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 4, false> >);
    static_assert(std::is_same_v< decltype(get<10>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, 4, false> >);
    static_assert(std::is_same_v< decltype(get<11>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 5, true> >);
    static_assert(std::is_same_v< decltype(get<12>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 5, true> >);
    static_assert(std::is_same_v< decltype(get<13>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitEdge, 6, false> >);
    static_assert(std::is_same_v< decltype(get<14>(e)), TypeDfsSearchEvent<DfsSearchEventType::VisitVertex, 5, false> >);
    static_assert(std::is_same_v< decltype(get<15>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, 5, false> >);
    static_assert(std::is_same_v< decltype(get<16>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 6, false> >);
    static_assert(std::is_same_v< decltype(get<17>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, 4, false> >);
    static_assert(std::is_same_v< decltype(get<18>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 4, false> >);
    static_assert(std::is_same_v< decltype(get<19>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, 3, false> >);
    static_assert(std::is_same_v< decltype(get<20>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 3, false> >);
    static_assert(std::is_same_v< decltype(get<21>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, 2, false> >);
    static_assert(std::is_same_v< decltype(get<22>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 1, false> >);
    static_assert(std::is_same_v< decltype(get<23>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, 1, false> >);
    static_assert(std::is_same_v< decltype(get<24>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackEdge, 0, false> >);
    static_assert(std::is_same_v< decltype(get<25>(e)), TypeDfsSearchEvent<DfsSearchEventType::BacktrackVertex, 0, true> >);
}

//
// DfsEdgeList
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueDfsEdgeList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    auto edges = ValueEdgeList{smarts};
    auto degrees = ValueVertexDegree{smarts, edges};
    auto incidentList = ValueIncidentList{smarts, edges, degrees};

    constexpr auto dfsEdges = ValueDfsEdgeList{smarts, incidentList}.data;

    static_assert(dfsEdges.size() == 4);
    static_assert(dfsEdges[0] == ValueDfsEdge{0, 0, 1, false});
    static_assert(dfsEdges[1] == ValueDfsEdge{1, 1, 2, false});
    static_assert(dfsEdges[2] == ValueDfsEdge{2, 1, 3, false});
    static_assert(dfsEdges[3] == ValueDfsEdge{3, 3, 0, true});
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeDfsEdgeList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    auto edges = TypeEdgeList{smarts};
    auto degrees = TypeVertexDegree{smarts, edges};
    auto incidentList = TypeIncidentList{smarts, edges, degrees};

    constexpr auto dfsEdges = TypeDfsEdgeList{smarts, incidentList}.data;

    static_assert(ctll::size(dfsEdges) == 4);
    static_assert(std::is_same_v< decltype(get<0>(dfsEdges)), TypeDfsEdge<0, 0, 1, false> >);
    static_assert(std::is_same_v< decltype(get<1>(dfsEdges)), TypeDfsEdge<1, 1, 2, false> >);
    static_assert(std::is_same_v< decltype(get<2>(dfsEdges)), TypeDfsEdge<2, 1, 3, false> >);
    static_assert(std::is_same_v< decltype(get<3>(dfsEdges)), TypeDfsEdge<3, 3, 0, true> >);
}

//
// Cyclemembership
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueCyclemembership")
{
    auto smarts = Smarts<"**1*(*)*1*">{};

    auto edges = ValueEdgeList{smarts};
    auto degrees = ValueVertexDegree{smarts, edges};
    auto incidentList = ValueIncidentList{smarts, edges, degrees};

    constexpr auto cycleMembership = ValueCycleMembership{smarts, incidentList}.data;
    constexpr auto v = std::get<0>(cycleMembership);
    constexpr auto e = std::get<1>(cycleMembership);

    static_assert(v[0] == false);
    static_assert(v[1] == true);
    static_assert(v[2] == true);
    static_assert(v[3] == false);
    static_assert(v[4] == true);
    static_assert(v[5] == false);

    static_assert(e[0] == false);
    static_assert(e[1] == true);
    static_assert(e[2] == false);
    static_assert(e[3] == true);
    static_assert(e[4] == true);
    static_assert(e[5] == false);
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeCyclemembership")
{
    auto smarts = Smarts<"**1*(*)*1*">{};

    auto edges = TypeEdgeList{smarts};
    auto degrees = TypeVertexDegree{smarts, edges};
    auto incidentList = TypeIncidentList{smarts, edges, degrees};

    constexpr auto cycleMembership = TypeCycleMembership{smarts, incidentList}.data;
    constexpr auto v = std::get<0>(cycleMembership);
    constexpr auto e = std::get<1>(cycleMembership);

    static_assert(v[0] == false);
    static_assert(v[1] == true);
    static_assert(v[2] == true);
    static_assert(v[3] == false);
    static_assert(v[4] == true);
    static_assert(v[5] == false);

    static_assert(e[0] == false);
    static_assert(e[1] == true);
    static_assert(e[2] == false);
    static_assert(e[3] == true);
    static_assert(e[4] == true);
    static_assert(e[5] == false);
}


//
// DfsBondList
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("ValueDfsBondList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    auto edges = ValueEdgeList{smarts};
    auto degrees = ValueVertexDegree{smarts, edges};
    auto incidentList = ValueIncidentList{smarts, edges, degrees};
    auto dfsEdges = ValueDfsEdgeList{smarts, incidentList};
    auto cycleMembership = ValueCycleMembership{smarts, incidentList};

    constexpr auto dfsBonds = ValueDfsBondList{smarts, dfsEdges, cycleMembership}.data;

    //identify_type<decltype(get<0>(dfsBonds))>{};

    using Atom0 = Atom<0, AnyAtom>;
    using Atom1 = Atom<1, AnyAtom>;
    using Atom2 = Atom<2, AnyAtom>;
    using Atom3 = Atom<3, AnyAtom>;

    static_assert(ctll::size(dfsBonds) == 4);
    static_assert(std::is_same_v< decltype(get<0>(dfsBonds)), DfsBond<Atom0, Atom1, ImplicitBond, true, false> >);
    static_assert(std::is_same_v< decltype(get<1>(dfsBonds)), DfsBond<Atom1, Atom2, ImplicitBond, false, false> >);
    static_assert(std::is_same_v< decltype(get<2>(dfsBonds)), DfsBond<Atom1, Atom3, ImplicitBond, true, false> >);
    static_assert(std::is_same_v< decltype(get<3>(dfsBonds)), DfsBond<Atom3, Atom0, ImplicitBond, true, true> >);
}

#endif // KITIMAR_VALUE_BASED

TEST_CASE("TypeDfsBondList")
{
    auto smarts = Smarts<"*1*(*)*1">{};

    auto edges = TypeEdgeList{smarts};
    auto degrees = TypeVertexDegree{smarts, edges};
    auto incidentList = TypeIncidentList{smarts, edges, degrees};
    auto dfsEdges = TypeDfsEdgeList{smarts, incidentList};
    auto cycleMembership = TypeCycleMembership{smarts, incidentList};

    constexpr auto dfsBonds = TypeDfsBondList{smarts, dfsEdges, cycleMembership}.data;

    //identify_type<decltype(get<0>(dfsBonds))>{};

    using Atom0 = Atom<0, AnyAtom>;
    using Atom1 = Atom<1, AnyAtom>;
    using Atom2 = Atom<2, AnyAtom>;
    using Atom3 = Atom<3, AnyAtom>;

    static_assert(ctll::size(dfsBonds) == 4);
    static_assert(std::is_same_v< decltype(get<0>(dfsBonds)), DfsBond<Atom0, Atom1, ImplicitBond, true, false> >);
    static_assert(std::is_same_v< decltype(get<1>(dfsBonds)), DfsBond<Atom1, Atom2, ImplicitBond, false, false> >);
    static_assert(std::is_same_v< decltype(get<2>(dfsBonds)), DfsBond<Atom1, Atom3, ImplicitBond, true, false> >);
    static_assert(std::is_same_v< decltype(get<3>(dfsBonds)), DfsBond<Atom3, Atom0, ImplicitBond, true, true> >);
}

//
// OptimizeIncidentList
//

#ifdef KITIMAR_VALUE_BASED

TEST_CASE("OptimizeIncidentList")
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
