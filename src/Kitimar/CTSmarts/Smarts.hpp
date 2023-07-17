#pragma once

#include "Algorithm/EdgeList.hpp"
#include "Algorithm/VertexDegree.hpp"
#include "Algorithm/IncidentList.hpp"
#include "Algorithm/DfsSearch.hpp"
#include "Algorithm/DfsSearchEvents.hpp"
#include "Algorithm/DfsEdgeList.hpp"
#include "Algorithm/CycleMembership.hpp"
#include "Algorithm/DfsBondList.hpp"

#include "Util.hpp"
#include "Parser/Actions.hpp"

#include <ctll/parser.hpp>

#include <ranges>
#include <algorithm>


namespace Kitimar::CTSmarts {

    template <ctll::fixed_string SMARTS, bool IgnoreInvalid = false>
    struct Smarts;



    constexpr auto cycleRank(auto numVertices, auto numEdges, auto numComponents)
    {
        return numEdges - numVertices + numComponents;
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
        auto vertexDegree = VertexDegree(smarts, edgeList);

        if constexpr (smarts.numAtoms < 3 || smarts.numAtoms != smarts.numBonds + 1)
            return -1;
        auto centralAtom = -1;
        for (std::size_t i = 0; i < smarts.numAtoms; ++i) {
            switch (vertexDegree.get(i)) {
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
