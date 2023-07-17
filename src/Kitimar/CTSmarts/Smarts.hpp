#pragma once

#include "Util.hpp"
#include "SmartsActions.hpp"

#include <ctll/parser.hpp>

#include <ranges>
#include <algorithm>

#include <iostream>

namespace Kitimar::CTSmarts {

    template <ctll::fixed_string SMARTS, bool IgnoreInvalid = false>
    struct Smarts;



    constexpr auto cycleRank(auto numVertices, auto numEdges, auto numComponents)
    {
        return numEdges - numVertices + numComponents;
    }


    //
    // EdgeList
    //

    struct Edge
    {
        int source = -1;
        int target = -1;
    };

    template<typename SmartsT>
    constexpr auto makeEdgeList(ctll::empty_list)
    {
        return std::array<Edge, SmartsT::numBonds>{};
    }

    template<typename SmartsT, typename Bond, typename ...Bonds>
    constexpr auto makeEdgeList(ctll::list<Bond, Bonds...> bonds)
    {
        auto [bond, tail] = ctll::pop_and_get_front(bonds);
        auto edges = makeEdgeList<SmartsT>(tail);
        auto bondIndex = SmartsT::numBonds - ctll::size(bonds);
        edges[bondIndex] = Edge{bond.source, bond.target};
        return edges;
    }

    template<typename SmartsT>
    struct EdgeList
    {
        // store source and target atom index for each edge
        static constexpr inline auto data = makeEdgeList<SmartsT>(SmartsT::bonds);

        static constexpr inline Edge get(int index) noexcept
        {
            return data[index];
        }

        constexpr EdgeList() noexcept {}
        constexpr EdgeList(SmartsT) noexcept {}
    };

    //
    // DegreeList
    //

    template<typename SmartsT, typename EdgeListT>
    constexpr auto makeDegreeList()
    {
        std::array<int, SmartsT::numAtoms> degrees = {};
        for (const auto &edge : EdgeListT::data) {
            ++degrees[edge.source];
            ++degrees[edge.target];
        }
        return degrees;
    }

    template<typename SmartsT, typename EdgeListT>
    struct DegreeList
    {
        // store degree for each vertex
        static constexpr inline auto data = makeDegreeList<SmartsT, EdgeListT>();

        static constexpr inline int get(int index) noexcept
        {
            return data[index];
        }

        static constexpr auto max() noexcept
        {
            return *std::ranges::max_element(data);
        }

        constexpr DegreeList() noexcept {}
        constexpr DegreeList(SmartsT, EdgeListT) noexcept {}
    };


    //
    // AdjacencyList
    //

    template<typename SmartsT, typename EdgeListT, typename DegreeListT>
    constexpr auto makeAdjacencyList()
    {
        constexpr auto stride = DegreeListT::max();
        std::array<int, SmartsT::numAtoms * stride> adj = {};

        for (auto &i : adj)
            i = -1; // FIXME: needed?

        std::array<int, SmartsT::numAtoms> sizes = {};
        for (auto i = 0; i < SmartsT::numBonds; ++i) {
            auto edge = EdgeListT::get(i);
            auto source = edge.source;
            auto target = edge.target;
            adj[stride * source + sizes[source]] = i;
            adj[stride * target + sizes[target]] = i;
            ++sizes[source];
            ++sizes[target];
        }

        return adj;
    }

    template<typename SmartsT, typename EdgeListT, typename DegreeListT>
    struct AdjacencyList
    {
        // store adjacent (or incident) bond indices for each vertex
        static constexpr inline auto data = makeAdjacencyList<SmartsT, EdgeListT, DegreeListT>();
        static constexpr inline auto edges = EdgeListT{};
        static constexpr inline auto degrees = DegreeListT{};
        static constexpr inline auto stride = DegreeListT::max();

        static constexpr auto get(int AtomIndex, int AdjIndex)
        {
            return data[stride * AtomIndex + AdjIndex];
        }

        constexpr AdjacencyList() noexcept {}
        constexpr AdjacencyList(SmartsT, EdgeListT, DegreeListT) noexcept {}
    };


    //
    // dfsSearch
    //

    template<typename SmartsT, typename Visitor, typename AdjacencyListT, int SourceIndex, int AdjIndex>
    constexpr void dfsSearchHelper(std::array<bool, SmartsT::numAtoms> &visitedVertices,
                                   std::array<bool, SmartsT::numBonds> &visitedEdges,
                                   Visitor &visitor)
    {
        constexpr auto sourceDegree = AdjacencyListT::degrees.get(SourceIndex);
        if constexpr (AdjIndex < sourceDegree) {
            constexpr auto edgeIndex = AdjacencyListT::get(SourceIndex, AdjIndex);

            if (visitedEdges[edgeIndex]) {
                dfsSearchHelper<SmartsT, Visitor, AdjacencyListT, SourceIndex, AdjIndex + 1>(visitedVertices, visitedEdges, visitor);
                return;
            }

            constexpr auto edge = AdjacencyListT::edges.get(edgeIndex);
            constexpr auto targetIndex = edge.source == SourceIndex ? edge.target : edge.source;
            auto isNewComponent = !visitedVertices[SourceIndex];
            auto isClosure = visitedVertices[targetIndex];

            visitor.visit(edgeIndex, SourceIndex, targetIndex, isNewComponent, isClosure);

            visitedEdges[edgeIndex] = true;
            visitedVertices[SourceIndex] = true;
            visitedVertices[targetIndex] = true;

            // dfs
            if (!isClosure)
                dfsSearchHelper<SmartsT, Visitor, AdjacencyListT, targetIndex, 0>(visitedVertices, visitedEdges, visitor);

            visitor.backtrack(edgeIndex, targetIndex, isClosure);

            // next incident bond
            dfsSearchHelper<SmartsT, Visitor, AdjacencyListT, SourceIndex, AdjIndex + 1>(visitedVertices, visitedEdges, visitor);
        } else if constexpr (SourceIndex == 0) {
            visitor.backtrack(SourceIndex);
        }
    }

    template<typename SmartsT, typename Visitor, typename AdjacencyListT>
    constexpr void dfsSearch(SmartsT, Visitor &visitor, AdjacencyListT)
    {
        std::array<bool, SmartsT::numAtoms> visitedAtoms = {};
        std::array<bool, SmartsT::numBonds> visitedBonds = {};
        dfsSearchHelper<SmartsT, Visitor, AdjacencyListT, 0, 0>(visitedAtoms, visitedBonds, visitor);
    }

    //
    // DfsSearchEvents
    //

    template<typename SmartsT>
    struct DfsSearchEventsVisitor
    {
        struct Event
        {
            enum Type {
                Invalid,
                VisitVertex, // index = vertex index, flag = isNewComponent
                VisitEdge, // index = edge index, flag = isClosure
                BacktrackVertex, // index = vertex index, flag = isEndComponent
                BacktrackEdge // index = edge index, flag = isClosure
            };

            Type type = Invalid;
            int index = -1;
            bool flag = false;
        };

        friend std::ostream& operator<<(std::ostream &os, const Event &event)
        {
            switch (event.type) {
                case Event::VisitVertex:
                    os << "VisitVertex( index = " << event.index << ", isNewComponent = " << event.flag << " )";
                    break;
                case Event::VisitEdge:
                    os << "VisitEdge( index = " << event.index << ", isClosure = " << event.flag << " )";
                    break;
                case Event::BacktrackVertex:
                    os << "BacktrackVertex( index = " << event.index << ", isEndComponent = " << event.flag << " )";
                    break;
                case Event::BacktrackEdge:
                    os << "BacktrackEdge( index = " << event.index << ", isClosure = " << event.flag << " )";
                    break;
                default:
                    os << "InvalidEvent";
                    break;
            }

            return os;
        }

        using Events = std::array<Event, 2 * (SmartsT::numAtoms + SmartsT::numBonds)>;

        constexpr DfsSearchEventsVisitor() noexcept {}
        constexpr DfsSearchEventsVisitor(SmartsT) noexcept {}

        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
        {
            if (isNewComponent)
                events[nextEventIndex++] = Event{Event::VisitVertex, source, true};
            events[nextEventIndex++] = Event{Event::VisitEdge, edge, isClosure};
            if (!isClosure)
                events[nextEventIndex++] = Event{Event::VisitVertex, target, false};
        }

        constexpr void backtrack(int edge, int target, bool isClosure) noexcept
        {
            if (!isClosure)
                events[nextEventIndex++] = Event{Event::BacktrackVertex, target, false};
            events[nextEventIndex++] = Event{Event::BacktrackEdge, edge, isClosure};
        }

        constexpr void backtrack(int source) noexcept
        {
            events[nextEventIndex++] = Event{Event::BacktrackVertex, source, true};
        }

        Events events = {};
        int nextEventIndex = 0;
    };

    template<typename SmartsT, typename AdjacencyListT>
    constexpr auto makeDfsSearchEvents()
    {
        DfsSearchEventsVisitor<SmartsT> visitor;
        dfsSearch(SmartsT{}, visitor, AdjacencyListT{});
        return visitor.events;
    }


    template<typename SmartsT, typename AdjacencyListT>
    struct DfsSearchEvents
    {
        static constexpr inline auto events = makeDfsSearchEvents<SmartsT, AdjacencyListT>();

        constexpr DfsSearchEvents() noexcept {}
        constexpr DfsSearchEvents(SmartsT, AdjacencyListT) noexcept {}
    };



    //
    // DfsEdgeList
    //

    struct DfsEdge : Edge
    {
        int index = -1;
        bool closure = false;
    };

    template<typename SmartsT>
    struct DfsEdgeListVisitor
    {

        constexpr DfsEdgeListVisitor() noexcept {}
        constexpr DfsEdgeListVisitor(SmartsT) noexcept {}

        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
        {
            edges[nextEdgeIndex++] = DfsEdge{{source, target}, edge, isClosure};
        }

        constexpr void backtrack(int edge, int target, bool isClosure) noexcept {}

        constexpr void backtrack(int source) noexcept {}


        std::array<DfsEdge, SmartsT::numBonds> edges;
        int nextEdgeIndex = 0;
    };

    template<typename SmartsT, typename AdjacencyListT>
    constexpr auto makeDfsEdgeList()
    {
        DfsEdgeListVisitor<SmartsT> visitor;
        dfsSearch(SmartsT{}, visitor, AdjacencyListT{});
        return visitor.edges;
    }


    template<typename SmartsT, typename AdjacencyListT>
    struct DfsEdgeList
    {
        static constexpr inline auto data = makeDfsEdgeList<SmartsT, AdjacencyListT>();

        static constexpr inline DfsEdge get(int index) noexcept
        {
            return data[index];
        }

        constexpr DfsEdgeList() noexcept {}
        constexpr DfsEdgeList(SmartsT, AdjacencyListT) noexcept {}
    };


    //
    // CycleMembership
    //

    template<typename SmartsT, typename AdjacencyListT>
    struct CycleMembershipVisitor
    {

        constexpr CycleMembershipVisitor() noexcept {}
        constexpr CycleMembershipVisitor(SmartsT) noexcept {}

        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
        {
            path[depth++] = edge;

            if (isClosure) {
                for (auto i = depth - 1; i >= 0; --i) {
                    auto e = AdjacencyListT::edges.get(path[i]);
                    edges[path[i]] = true;
                    vertices[e.source] = true;
                    vertices[e.target] = true;
                    if (i < depth - 1)
                        if (e.source == target || e.target == target)
                            break;
                }
            }
        }

        constexpr void backtrack(int edge, int target, bool isClosure) noexcept
        {
            --depth;
        }

        constexpr void backtrack(int source) noexcept {}

        std::array<int, SmartsT::numBonds> path = {};
        std::array<bool, SmartsT::numAtoms> vertices = {};
        std::array<bool, SmartsT::numBonds> edges = {};
        int depth = 0;
    };

    template<typename SmartsT, typename AdjacencyListT>
    constexpr auto makeCycleMembership()
    {
        CycleMembershipVisitor<SmartsT, AdjacencyListT> visitor;
        dfsSearch(SmartsT{}, visitor, AdjacencyListT{});
        return std::make_tuple(visitor.vertices, visitor.edges);
    }

    template<typename SmartsT, typename EdgeListT>
    struct CycleMembership
    {
        static constexpr inline auto data = makeCycleMembership<SmartsT, EdgeListT>();
        static constexpr inline auto vertices = std::get<0>(data);
        static constexpr inline auto edges = std::get<1>(data);

        constexpr CycleMembership() noexcept {}
        constexpr CycleMembership(SmartsT, EdgeListT) noexcept {}
    };



    //
    // DfsBondList
    //

    template<int Source, int Target, bool IsCyclic, bool IsRingClosure, typename SourceExpr, typename TargetExpr, typename BondExpr>
    struct DfsBond
    {
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
        static constexpr inline auto sourceExpr = SourceExpr();
        static constexpr inline auto targetExpr = TargetExpr();
        static constexpr inline auto bondExpr = BondExpr();
    };

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT, int DfsEdgeIndex>
    constexpr auto makeDfsBondListHelper()
    {
        if constexpr (DfsEdgeIndex == SmartsT::numBonds) {
            return ctll::empty_list{};
        } else {
            constexpr auto edge = DfsEdgeListT::get(DfsEdgeIndex);

            constexpr auto bondExpr = get<edge.index>(SmartsT::bonds).expr;
            auto sourceExpr = get<edge.source>(SmartsT::atoms);
            auto targetExpr = get<edge.target>(SmartsT::atoms);
            constexpr auto isCyclic = CycleMembershipT::edges[edge.index];

            //auto dfsBond = DfsBond<edge.source, edge.target, false, edge.closure, decltype(sourceExpr), decltype(targetExpr), decltype(bondExpr)>{};
            auto dfsBond = DfsBond<edge.source, edge.target, isCyclic, edge.closure, decltype(sourceExpr), decltype(targetExpr), decltype(bondExpr)>{};

            return ctll::push_front(dfsBond, makeDfsBondListHelper<SmartsT, DfsEdgeListT, CycleMembershipT, DfsEdgeIndex + 1>());
        }
    }

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    constexpr auto makeDfsBondList()
    {
        return makeDfsBondListHelper<SmartsT, DfsEdgeListT, CycleMembershipT, 0>();
    }


    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    struct DfsBondList
    {
        static constexpr inline auto data = makeDfsBondList<SmartsT, DfsEdgeListT, CycleMembershipT>();

        static constexpr inline DfsEdge get(int index) noexcept
        {
            return data[index];
        }

        constexpr DfsBondList() noexcept {}
        constexpr DfsBondList(SmartsT, DfsEdgeListT, CycleMembershipT) noexcept {}
    };













    /*

    struct DfsSearch
    {

        static constexpr inline auto NoAtomVisitor = [] (auto, auto, auto ctx) { return ctx; };
        static constexpr inline auto NoBondVisitor = [] (auto, auto, auto, auto, auto, auto ctx) { return ctx; };
        static constexpr inline auto NoAtomBacktrack = [] (auto ctx) { return ctx; };
        static constexpr inline auto NoBondBacktrack = [] (auto ctx) { return ctx; };

        static constexpr auto visitAtom(auto smarts, auto atomIdx, auto atomVisitor, auto ctx, auto visitedAtoms) noexcept
        {
            if constexpr (!ctll::exists_in(atomIdx, visitedAtoms))
                return atomVisitor(smarts, atomIdx, ctx);
            else
                return ctx;
        }

        static constexpr auto backtrackAtom(auto cond, auto atomBacktrack, auto ctx) noexcept
        {
            if constexpr (cond.value)
                return atomBacktrack(ctx);
            else
                return ctx;
        }

        template<int sourceIdx, int adjIdx, typename AtomVisitor, typename BondVisitor, typename Context, typename VB, typename VA>
        static constexpr auto visit(auto smarts, auto adjList, AtomVisitor atomVisitor, BondVisitor bondVisitor, auto atomBacktrack, auto bondBacktrack, Context ctx, VB visitedBonds, VA visitedAtoms) noexcept
        {
            //auto adjBondIdxs = get<sourceIdx>(smarts.adjList);
            if constexpr (adjIdx >= adjList.degrees.get(sourceIdx)) {
                // A leaf node has been reached
                return std::make_tuple(visitedBonds, visitedAtoms, ctx);
            } else {
                //auto bondIdx = get<adjIdx>(adjBondIdxs);
                constexpr auto bondIdx = std::integral_constant<int, adjList.get(sourceIdx, adjIdx)>{};
                if constexpr (ctll::exists_in(bondIdx, visitedBonds)) {
                    // Skip visited bonds
                    return visit<sourceIdx, adjIdx + 1>(smarts, adjList, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx, visitedBonds, visitedAtoms);;
                } else {
                    auto bond = get<bondIdx.value>(smarts.bonds);
                    constexpr auto targetIdx = bond.source == sourceIdx ? bond.target : bond.source;
                    //auto targetExpr = get<targetIdx>(smarts.atoms);
                    constexpr auto isRingClosure = ctll::exists_in(std::integral_constant<int, targetIdx>(), visitedAtoms);
                    // Visit atoms and bond
                    auto ctx2 = visitAtom(smarts, std::integral_constant<int, sourceIdx>{}, atomVisitor, ctx, visitedAtoms);
                    auto ctx3  = bondVisitor(smarts, std::integral_constant<int, sourceIdx>{}, std::integral_constant<int, targetIdx>{}, bond.expr, std::integral_constant<bool, isRingClosure>{}, ctx2);
                    auto ctx4 = visitAtom(smarts, std::integral_constant<int, targetIdx>{}, atomVisitor, ctx3, visitedAtoms);
                    // Mark bond and atoms as visited
                    auto visitedBonds2 = ctll::push_front(bondIdx, visitedBonds);
                    auto visitedAtoms2 = ctll::push_front(std::integral_constant<int, sourceIdx>{}, ctll::push_front(std::integral_constant<int, targetIdx>{}, visitedAtoms));
                    // Recursive search....
                    auto [visitedBonds3, visitedAtoms3, ctx5] = visit<targetIdx, 0>(smarts, adjList, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx4, visitedBonds2, visitedAtoms2);
                    auto ctx6 = backtrackAtom(std::integral_constant<bool, !isRingClosure>{}, atomBacktrack, ctx5);
                    auto ctx7 = bondBacktrack(ctx6);
                    // bond to next nbr
                    return visit<sourceIdx, adjIdx + 1>(smarts, adjList, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx7, visitedBonds3, visitedAtoms3);
                }
            }
        }

        template<typename Context = ctll::empty_list>
        static constexpr auto visit(auto smarts, auto adjList, auto atomVisitor, auto bondVisitor, auto atomBacktrack, auto bondBacktrack, Context ctx = {}) noexcept
        {
            auto ctx2 = std::get<2>(visit<0, 0>(smarts, adjList, atomVisitor, bondVisitor, atomBacktrack, bondBacktrack, ctx, ctll::empty_list{}, ctll::empty_list{}));
            return atomBacktrack(ctx2);
        }

        template<typename Context = ctll::empty_list>
        static constexpr auto visit(auto smarts, auto adjList, auto atomVisitor, auto bondVisitor, Context ctx = {}) noexcept
        {
            return visit(smarts, adjList, atomVisitor, bondVisitor, NoAtomBacktrack, NoBondBacktrack, ctx);
        }

    };

    */

    /*
    constexpr auto getDfsAtoms(auto smarts, auto adjList) noexcept
    {
        constexpr auto atomVisitor = [] (auto smarts, auto atomIdx, auto ctx) {
            return ctll::push_front(atomIdx, ctx);
        };
        return ctll::rotate(DfsSearch::visit(smarts, adjList, atomVisitor, DfsSearch::NoBondVisitor));
    }
    */

    //
    // Cycle-membership
    //

    /*
    template<int BondIdx, typename Atoms, typename Bonds, typename Path>
    struct CycleMembershipContext
    {
        static constexpr inline auto atoms = Atoms{};
        static constexpr inline auto bonds = Bonds{};
        static constexpr inline auto path = Path{};
        static constexpr inline auto bondIdx = BondIdx;

        constexpr CycleMembershipContext() noexcept = default;
        constexpr CycleMembershipContext(Atoms, Bonds, Path, Number<BondIdx>) noexcept {}

    };

    template<int SourceIdx, int TargetIdx, int BondIdx>
    struct CylceMembershipBond
    {
        static constexpr inline auto source = SourceIdx;
        static constexpr inline auto target = TargetIdx;
        static constexpr inline auto idx = BondIdx;
    };


    constexpr auto getCycleMembershipHelper(auto targetIdx, ctll::empty_list path, auto ctx) noexcept
    {
        return ctx;
    }

    template<typename Bond, typename ...Bonds>
    constexpr auto getCycleMembershipHelper(auto targetIdx, ctll::list<Bond, Bonds...> path, auto ctx) noexcept
    {
        auto atoms = ctll::add_item(Number<Bond::target>{}, ctll::add_item(Number<Bond::source>{}, ctx.atoms));
        auto bonds = ctll::add_item(Number<Bond::idx>{}, ctx.bonds);
        auto ctx2 = CycleMembershipContext{atoms, bonds, ctx.path, Number<ctx.bondIdx>{}};
        if constexpr (Bond::idx != ctx.bondIdx && (Bond::source == targetIdx.value || Bond::target == targetIdx.value))
            return ctx2;
        else
            return getCycleMembershipHelper(targetIdx, ctll::list<Bonds...>{}, ctx2);
    }

    constexpr auto getCycleMembership(auto smarts, auto adjList) noexcept
    {
        constexpr auto bondVisitor = [] (auto smarts, auto sourceIdx, auto targetIdx, auto expr, auto isRingClosure, auto ctx) {
            auto path = ctll::push_front(CylceMembershipBond<sourceIdx.value, targetIdx.value, ctx.bondIdx>{}, ctx.path);

            if constexpr (isRingClosure.value) {
                auto ctx2 = getCycleMembershipHelper(targetIdx, path, ctx);
                return CycleMembershipContext(ctx2.atoms, ctx2.bonds, path, Number<ctx2.bondIdx+1>{});
            } else {
                return CycleMembershipContext(ctx.atoms, ctx.bonds, path, Number<ctx.bondIdx+1>{});
            }
        };

        constexpr auto bondBacktrack = [] (auto ctx) {
            auto path = ctll::pop_front(ctx.path);
            return CycleMembershipContext{ctx.atoms, ctx.bonds, path, Number<ctx.bondIdx>{}};
        };

        auto ctx = CycleMembershipContext<0, ctll::empty_list, ctll::empty_list, ctll::empty_list>{};
        return DfsSearch::visit(smarts, adjList, DfsSearch::NoAtomVisitor, bondVisitor, DfsSearch::NoAtomBacktrack, bondBacktrack, ctx);
    }

    template<int BondIdx>
    constexpr auto isCyclic(ctll::empty_list)
    {
        return false;
    }

    template<int BondIdx, typename Bond, typename ...Bonds>
    constexpr auto isCyclic(ctll::list<Bond, Bonds...>)
    {
        if (Bond::value == BondIdx)
            return true;
        return isCyclic<BondIdx>(ctll::list<Bonds...>{});
    }

    template<int BondIdx>
    constexpr auto addCycleMembership(ctll::empty_list, auto cyclicBondIdxs) noexcept
    {
        return ctll::empty_list{};
    }

    template<int BondIdx, int S, int T, bool IC, bool IRC, typename SE, typename TE, typename BE, typename ...CycleBonds>
    constexpr auto addCycleMembership(ctll::list<DfsBond<S, T, IC, IRC, SE, TE, BE>, CycleBonds...>, auto cyclicBondIdxs) noexcept
    {
        auto dfsBond = DfsBond<S, T, isCyclic<BondIdx>(cyclicBondIdxs), IRC, SE, TE, BE>{};
        return ctll::push_front(dfsBond, addCycleMembership<BondIdx+1>(ctll::list<CycleBonds...>{}, cyclicBondIdxs));
    }
    */

    /*
    constexpr auto getDfsBonds(auto smarts, auto adjList) noexcept
    {
        constexpr auto bondVisitor = [] (auto smarts, auto sourceIdx, auto targetIdx, auto expr, auto isRingClosure, auto ctx) {
            auto sourceExpr = get<sourceIdx.value>(smarts.atoms);
            auto targetExpr = get<targetIdx.value>(smarts.atoms);
            return ctll::push_front(DfsBond<sourceIdx.value, targetIdx.value, false, isRingClosure.value,
                                            decltype(sourceExpr), decltype(targetExpr), decltype(expr)>(), ctx);
        };
        auto dfsBonds = ctll::rotate(DfsSearch::visit(smarts, adjList, DfsSearch::NoAtomVisitor, bondVisitor));
        return dfsBonds;
        //return addCycleMembership<0>(dfsBonds, getCycleMembership(smarts, adjList).bonds);
    }
    */







    /*
    template<typename SmartsT, typename Bonds = decltype(SmartsT::bonds)>
    constexpr auto getDfsBonds(SmartsT smarts, Bonds bonds = {}) noexcept
    {
        if constexpr (ctll::empty(bonds)) {
            return ctll::empty_list{};
        } else {
            auto [bond, tail] = ctll::pop_and_get_front(bonds);
            auto bondIndex = smarts.numBonds - ctll::size(bonds);
            auto sourceExpr = get<bond.source>(smarts.atoms);
            auto targetExpr = get<bond.target>(smarts.atoms);

            auto dfsBond = DfsBond<bond.source, bond.target, false, false, decltype(sourceExpr), decltype(targetExpr), decltype(bond.expr)>();

            return ctll::push_front(dfsBond, getDfsBonds(smarts, tail));
        }
    }
    */




    //
    // Captures
    //

    template<int N, typename Classes>
    constexpr auto captureIndex(Classes classes)
    {
        if constexpr (ctll::empty(classes)) {
            return -1;
        } else {
            auto cls = ctll::front(classes);
            if (cls.n == N)
                return cls.atomIndex;
            return captureIndex<N>(ctll::pop_front(classes));
        }
    }


    template<int I, typename Classes, typename Mapping>
    constexpr auto captureMappingHelper(Classes classes, Mapping mapping)
    {
        constexpr auto index = captureIndex<I>(classes);
        if constexpr (index < 0) {
            return ctll::rotate(mapping);
        } else {
            return captureMappingHelper<I + 1>(classes, ctll::push_front(Number<index>(), mapping));
        }
    }

    constexpr auto captureMapping(auto smarts)
    {
        constexpr auto mapping = captureMappingHelper<1>(smarts.context.params.classes, ctll::empty_list());
        return toArray(mapping);
    }

    //
    // Optimized cases
    //



    constexpr auto getCentralAtom(auto smarts)
    {
        auto edgeList = EdgeList(smarts); // FIXME
        auto degreeList = DegreeList(smarts, edgeList);

        if constexpr (smarts.numAtoms < 3 || smarts.numAtoms != smarts.numBonds + 1)
            return -1;
        auto centralAtom = -1;
        for (std::size_t i = 0; i < smarts.numAtoms; ++i) {
            switch (degreeList.get(i)) {
                case 0:
                    return -1;
                case 1:
                    continue;
                default:
                    if (centralAtom != -1)
                        return -1;
                    centralAtom = i;
                    break;
            }
        }
        return centralAtom;
    }




    //
    // Errors
    //

    template<int I>
    struct Pos {};

    struct NoError {};

    template<typename>
    struct EmptyBracketAtomError {};

    //template<typename>
    struct ConflicingRingBondError {};

    template<typename>
    constexpr auto UnmatchedRingBondError() { return false; }


    template<typename ...Ns>
    struct RingBondIds {};

    template<typename ...Ids>
    constexpr auto ringBondIds(ctll::empty_list, ctll::list<Ids...>)
    {
        return RingBondIds<Ids...>();
    }

    template<typename Bond, typename ...Bonds, typename Ids = ctll::empty_list>
    constexpr auto ringBondIds(ctll::list<Bond, Bonds...>, Ids = {})
    {
        return ringBondIds(ctll::list<Bonds...>(), ctll::push_front(Number<Bond::n>(), Ids()));
    }


    template<typename ...Expr>
    constexpr auto handleTotalH(ctll::list<Expr...> atoms)
    {
        return transform(atoms, [] (auto expr) {
            if constexpr (SmartsActions::isTotalHExpr(expr))
                return SmartsActions::toTotalHExpr(expr);
            else
                return expr;
        });
    }


    //template<typename Atoms, typename Bonds>
    template <ctll::fixed_string SMARTS, bool IgnoreInvalid>
    struct Smarts
    {
        using Result = ctll::parser<SmartsGrammar, SMARTS, SmartsActions>::template output<SmartsContext<>>;
        using Context = Result::output_type;

        static constexpr inline auto smarts = SMARTS;

        static constexpr auto input()
        {
            auto str = SMARTS | std::views::transform([] (auto c) { return static_cast<char>(c); });
            return std::string(str.begin(), str.end());
        }

        static constexpr inline auto context = Context();
        static constexpr inline auto valid = Result::is_correct;
        static constexpr inline auto position = Result::position;


        static constexpr inline auto atoms = handleTotalH(ctll::rotate(Context::atoms));
        static constexpr inline auto bonds = ctll::rotate(Context::bonds);
        static constexpr inline auto numAtoms = ctll::size(atoms);
        static constexpr inline auto numBonds = ctll::size(bonds);

        static constexpr inline auto isSingleAtom = numAtoms == 1;
        static constexpr inline auto isSingleBond = numAtoms == 2 && numBonds == 1;




        template<typename ErrorTag>
        static constexpr auto getError(ErrorTag)
        {
            if constexpr (std::is_same_v<ErrorTag, EmptyBracketAtomTag>)
                return EmptyBracketAtomError<Pos<position>>();
            else if constexpr (std::is_same_v<ErrorTag, ConflicingRingBondTag>)
                return ConflicingRingBondError{};
            else if constexpr (!ctll::empty(context.params.ringBonds))
                return UnmatchedRingBondError<decltype(ringBondIds(context.params.ringBonds))>();
            else
                return NoError();
        }

        static constexpr inline auto error = getError(context.params.error);





        //static_assert(IgnoreInvalid || valid);
        static_assert(IgnoreInvalid || std::is_same_v<const NoError, decltype(error)>); // FIXME
        //static_assert(IgnoreInvalid || error); // FIXME

    };

















} // namespace ctsmarts
