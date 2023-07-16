#pragma once

#include "Util.hpp"
#include "SmartsActions.hpp"

#include <ctll/parser.hpp>

#include <ranges>
#include <algorithm>

namespace Kitimar::CTSmarts {

    template <ctll::fixed_string SMARTS, bool IgnoreInvalid = false>
    struct Smarts;






    // EdgeList

    struct Edge
    {
        int source;
        int target;
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

        constexpr EdgeList(SmartsT) noexcept {}
    };

    // DegreeList

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

        constexpr DegreeList(SmartsT, EdgeListT) noexcept {}
    };


    // AdjacencyList


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
        static constexpr inline auto degrees = DegreeListT::data;
        static constexpr inline auto stride = DegreeListT::max();

        template<int AtomIndex, int AdjIndex>
        constexpr auto get(Number<AtomIndex>, Number<AdjIndex>)
        {
            return data[stride * AtomIndex + AdjIndex];
        }

        constexpr AdjacencyList(SmartsT, EdgeListT, DegreeListT) noexcept {}
    };













    //
    // Depth-first search bonds
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
            if constexpr (adjIdx >= adjList.degrees[sourceIdx]) {
                // A leaf node has been reached
                return std::make_tuple(visitedBonds, visitedAtoms, ctx);
            } else {
                //auto bondIdx = get<adjIdx>(adjBondIdxs);
                constexpr auto bondIdx = std::integral_constant<int, adjList.get(Number<sourceIdx>{}, Number<adjIdx>{})>{};
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

    constexpr auto getDfsAtoms(auto smarts, auto adjList) noexcept
    {
        constexpr auto atomVisitor = [] (auto smarts, auto atomIdx, auto ctx) {
            return ctll::push_front(atomIdx, ctx);
        };
        return ctll::rotate(DfsSearch::visit(smarts, adjList, atomVisitor, DfsSearch::NoBondVisitor));
    }

    //
    // Cycle-membership
    //

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

    constexpr auto getDfsBonds(auto smarts, auto adjList) noexcept
    {
        constexpr auto bondVisitor = [] (auto smarts, auto sourceIdx, auto targetIdx, auto expr, auto isRingClosure, auto ctx) {
            auto sourceExpr = get<sourceIdx.value>(smarts.atoms);
            auto targetExpr = get<targetIdx.value>(smarts.atoms);
            return ctll::push_front(DfsBond<sourceIdx.value, targetIdx.value, false, isRingClosure.value,
                                            decltype(sourceExpr), decltype(targetExpr), decltype(expr)>(), ctx);
        };
        auto dfsBonds = ctll::rotate(DfsSearch::visit(smarts, adjList, DfsSearch::NoAtomVisitor, bondVisitor));
        return addCycleMembership<0>(dfsBonds, getCycleMembership(smarts, adjList).bonds);
    }




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
