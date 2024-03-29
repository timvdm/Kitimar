#pragma once

#include <concepts>
#include <ranges>

namespace Kitimar::Molecule {

    template<typename R, typename AtomBond>
    concept AtomBondRange = std::ranges::input_range<R> &&
                        std::convertible_to<std::remove_cvref_t<std::ranges::range_value_t<R>>, AtomBond>;

    // FIXME atom != atom

    template<typename Mol>
    concept AtomList = requires (Mol &mol)
    {
        { num_atoms(mol) } -> std::convertible_to<std::size_t>;
        { get_atom(mol, 0) };
        { get_atoms(mol) } -> AtomBondRange<std::remove_cvref_t<decltype(get_atom(mol, 0))>>;
        { get_index(mol, get_atom(mol, 0)) } -> std::convertible_to<std::size_t>;
        { null_atom(mol) };
    };

    template<typename Mol>
    concept BondList = requires (Mol &mol)
    {
        { num_bonds(mol) } -> std::convertible_to<std::size_t>;
        { get_bond(mol, 0) };
        { get_bonds(mol) } -> AtomBondRange<std::remove_cvref_t<decltype(get_bond(mol, 0))>>;
        { get_index(mol, get_bond(mol, 0)) } -> std::convertible_to<std::size_t>;
        { get_source(mol, get_bond(mol, 0)) } -> std::convertible_to<decltype(get_atom(mol, 0))>;
        { get_target(mol, get_bond(mol, 0)) } -> std::convertible_to<decltype(get_atom(mol, 0))>;
        { null_bond(mol) };
    };

    template<typename Mol>
    concept MoleculeGraph = AtomList<Mol> &&
                            BondList<Mol>;

    template<typename Mol>
    concept IncidentBondList = requires (Mol &mol)
    {
        { get_degree(mol, get_atom(mol, 0)) } -> std::integral;
        { get_bonds(mol, get_atom(mol, 0)) } -> AtomBondRange<std::remove_cvref_t<decltype(get_bond(mol, 0))>>;
    };

    template<typename Mol>
    concept AdjacentAtomList = requires (Mol &mol)
    {
        { get_degree(mol, get_atom(mol, 0)) } -> std::integral;
        { get_nbrs(mol, get_atom(mol, 0)) } -> AtomBondRange<std::remove_cvref_t<decltype(get_atom(mol, 0))>>;
    };

    template<typename Mol>
    concept ElementLayer = requires (Mol &mol)
    {
        { get_element(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept IsotopeLayer = requires (Mol &mol)
    {
        { get_isotope(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept ChargeLayer = requires (Mol &mol)
    {
        { get_charge(mol, get_atom(mol, 0)) } -> std::signed_integral;
    };

    // bond order sum (including implicit hydrogens)
    template<typename Mol>
    concept ValenceLayer = requires (Mol &mol)
    {
        { get_valence(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept BondOrderLayer = requires (Mol &mol)
    {
        { get_order(mol, get_bond(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept ImplicitHydrogensLayer = requires (Mol &mol)
    {
        { get_implicit_hydrogens(mol, get_atom(mol, 0)) } -> std::integral;
        { get_total_hydrogens(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept MobileHydrogensLayer = requires (Mol &mol)
    {
        { get_mobile_hydrogens(mol, get_atom(mol, 0)) } -> std::integral;
    };

    template<typename Mol>
    concept RingLayer = requires (Mol &mol)
    {
        { is_ring_atom(mol, get_atom(mol, 0)) } -> std::same_as<bool>;
        { is_ring_bond(mol, get_bond(mol, 0)) } -> std::same_as<bool>;
    };

    template<typename Mol>
    concept RingSetLayer = requires (Mol &mol)
    {
        { is_in_ring_size(mol, get_atom(mol, 0), 6) } -> std::same_as<bool>;
        { get_ring_count(mol, get_atom(mol, 0), 6) } -> std::integral;
        { get_ring_degree(mol, get_atom(mol, 0), 6) } -> std::integral;

    };

    template<typename Mol>
    concept AromaticLayer = requires (Mol &mol)
    {
        { is_aromatic_atom(mol, get_atom(mol, 0)) } -> std::same_as<bool>;
        { is_aromatic_bond(mol, get_bond(mol, 0)) } -> std::same_as<bool>;
    };

    template<typename Mol>
    concept Molecule = MoleculeGraph<Mol> &&
                       IncidentBondList<Mol> &&
                       AdjacentAtomList<Mol> &&
                       ElementLayer<Mol> &&
                       IsotopeLayer<Mol> &&
                       ChargeLayer<Mol> &&
                       BondOrderLayer<Mol> &&
                       //ImplicitHydrogensLayer<Mol> &&
                       AromaticLayer<Mol>;

    template<Molecule Mol>
    struct MoleculeTraits
    {
        using Atom = decltype(get_atom(std::declval<Mol>(), 0));
        using Bond = decltype(get_bond(std::declval<Mol>(), 0));
        constexpr static inline bool SameAtomBondType = std::is_same_v<Atom, Bond>;
    };

    constexpr auto get_nbr(const auto &mol, const auto &bond, const auto &atom) noexcept -> decltype(get_atom(mol, 0))
    {
        auto source = get_source(mol, bond);
        return source != atom ? source : get_target(mol, bond);
    }

    constexpr auto get_bond(const auto &mol, auto source, auto target) noexcept -> decltype(get_bond(mol, 0))
    {
        auto sourceIndex1 = get_index(mol, source);
        auto targetIndex1 = get_index(mol, target);
        for (auto bond : get_bonds(mol, source)) {
            auto sourceIndex2 = get_index(mol, get_source(mol, bond));
            auto targetIndex2 = get_index(mol, get_target(mol, bond));
            if ((sourceIndex1 == sourceIndex2 && targetIndex1 == targetIndex2) ||
                (sourceIndex1 == targetIndex2 && targetIndex1 == sourceIndex2))
                return bond;
        }
        return null_bond(mol);
    }

} // namespace Kitimar::Molecule

#include <vector>
#include <algorithm>
#include <functional>
#include <type_traits>

namespace Kitimar::CTSmarts {

    namespace impl {

        /**
         * @brief Convert atom index map to atom map using capture set.
         */
        auto toCapture(const auto &mol, auto smarts,
                       const auto &captureSet, bool found, const auto &map)
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto NumCaptures = std::tuple_size_v<std::remove_cvref_t<decltype(captureSet)>>;
            if constexpr (NumCaptures) {
                std::array<Atom, NumCaptures> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0UL; i < NumCaptures; ++i)
                        atoms[i] = get_atom(mol, map[captureSet[i]]);
                return atoms;
            } else {
                std::array<Atom, smarts.numAtoms> atoms = {};
                if (!found)
                    atoms.fill(null_atom(mol));
                else
                    for (auto i = 0UL; i < smarts.numAtoms; ++i)
                        atoms[i] = get_atom(mol, map[i]); // FIXME: null atoms...
                return atoms;
            }
        }

        auto toCaptures(const auto &mol, const auto &iso,
                        const auto &captureSet, const auto &maps) noexcept
        {
            using Atom = decltype(get_atom(mol, 0));
            static constexpr auto NumCaptures = std::tuple_size_v<std::remove_cvref_t<decltype(captureSet)>>;
            static constexpr auto M = NumCaptures ? NumCaptures : decltype(iso.smarts)::numAtoms;
            auto r = maps | std::views::transform([&] (const auto &map) {
                return toCapture(mol, iso.smarts, captureSet, true, map);
            });
            return std::vector<std::array<Atom, M>>{r.begin(), r.end()};
        }

        auto captureMatchAtoms(const auto &mol, auto smarts,
                               const auto &captureSet, bool found, const auto &map)
        {
            return std::tuple_cat(std::make_tuple(found), toCapture(mol, smarts, captureSet, found, map));
        }

        constexpr std::size_t captureHash(const auto &mol, const auto &capture)
        {
            auto atoms = std::vector<bool>(num_atoms(mol));
            for (const auto &atom : capture)
                atoms[get_index(mol, atom)] = true;
            return std::hash<std::vector<bool>>()(atoms);
        }

        template<typename Proj = std::identity>
        constexpr auto numInversions(const auto &permutation, Proj proj = {})
        {
            std::size_t n = 0;
            for (auto i = 0; i < permutation.size(); ++i)
                for (auto j = i + 1; j < permutation.size(); ++j)
                    if (std::invoke(proj, permutation[i]) > std::invoke(proj, permutation[j]))
                        ++n;
            return n;
        }

        constexpr bool singleAtomMatch(auto smarts, const auto &mol, const auto &atom)
        {
            return matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr);
        }

        constexpr bool singleBondMatchHelper(auto smarts, const auto &mol, const auto &bond, const auto &source, const auto &target)
        {
            return matchAtomExpr(mol, source, get<0>(smarts.atoms).expr) && matchAtomExpr(mol, target, get<1>(smarts.atoms).expr);
        }

        // 0 -> no match
        // 1 -> source is SMARTS atom 0, target is SMARTS atom 1
        // 2 -> source is SMARTS atom 1, target is SMARTS atom 0
        constexpr int singleBondMatch(auto smarts, const auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return 1;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return 2;
            return 0;
        }

        constexpr int singleBondCount(auto smarts, const auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return 0;
            auto n = 0;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                ++n;
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                ++n;
            return n;
        }

        template<typename Map>
        constexpr auto singleBondMap(auto smarts, const auto &mol, const auto &bond)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return std::make_tuple(false, Map{});

            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                return std::make_tuple(true, Map{get_index(mol, get_source(mol, bond)),
                                                 get_index(mol, get_target(mol, bond))});
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                return std::make_tuple(true, Map{get_index(mol, get_target(mol, bond)),
                                                 get_index(mol, get_source(mol, bond))});

            return std::make_tuple(false, Map{});
        }

        template<typename Map>
        constexpr void singleBondMaps(auto smarts, const auto &mol, const auto &bond, std::vector<Map> &maps)
        {
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return;
            if (singleBondMatchHelper(smarts, mol, bond, source, target))
                maps.push_back({get_index(mol, get_source(mol, bond)),
                                get_index(mol, get_target(mol, bond))});
            if (singleBondMatchHelper(smarts, mol, bond, target, source))
                maps.push_back({get_index(mol, get_target(mol, bond)),
                                get_index(mol, get_source(mol, bond))});
        }

        constexpr auto singleBondCapture(auto smarts, const auto &mol, const auto &bond, int singleBondMatchType)
        {
            constexpr auto cap = captureMapping(smarts);
            if constexpr (cap.size() == 1) {
                switch (singleBondMatchType) {
                    case 1:
                        return cap[0] == 0 ? std::make_tuple(true, get_source(mol, bond)) :
                                             std::make_tuple(true, get_target(mol, bond));
                    case 2:
                        return cap[0] == 0 ? std::make_tuple(true, get_target(mol, bond)) :
                                             std::make_tuple(true, get_source(mol, bond));
                    default:
                        return std::make_tuple(false, null_atom(mol));
                }
            } else {
                if (!singleBondMatchType)
                    return std::make_tuple(false, null_atom(mol), null_atom(mol));

                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                if constexpr (cap.size())
                    if (cap[0] > cap[1])
                        std::swap(source, target);

                if (singleBondMatchType == 1)
                    return std::make_tuple(true, source, target);
                return std::make_tuple(true, target, source);
            }
        }

        template<typename AtomMap>
        constexpr void singleBondCaptures(auto smarts, const auto &mol, const auto &bond, const auto &captureSet, std::vector<AtomMap> &maps, bool unique)
        {
            using IndexMap = std::array<decltype(get_index(mol, get_atom(mol, 0))), 2>;
            auto source = get_source(mol, bond);
            auto target = get_target(mol, bond);
            if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                return;

            if (singleBondMatchHelper(smarts, mol, bond, source, target)) {
                auto map = IndexMap{get_index(mol, get_source(mol, bond)),
                                    get_index(mol, get_target(mol, bond))};
                maps.push_back(impl::toCapture(mol, smarts, captureSet, true, map));
                if (unique)
                    return;
            }

            if (singleBondMatchHelper(smarts, mol, bond, target, source)) {
                auto map = IndexMap{get_index(mol, get_target(mol, bond)),
                                    get_index(mol, get_source(mol, bond))};
                maps.push_back(impl::toCapture(mol, smarts, captureSet, true, map));
            }
        }

        constexpr auto centralAtomMap(auto smarts, const auto &mol, const auto &atom)
        {
            std::array<int, smarts.numAtoms> map;
        }

    } // namespace impl

// function_search(mol) -> function(mol, Search)

#define CTSMARTS_API_SEARCH(function, Search, search) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig> \
    constexpr auto function##_##search(const Molecule::Molecule auto &mol) \
    { return function<SMARTS, Config, SearchType::Search>(mol); }

#define CTSMARTS_API_UNIQUE(function) CTSMARTS_API_SEARCH(function, Unique, unique)
#define CTSMARTS_API_ALL(function) CTSMARTS_API_SEARCH(function, All, all)

// function_arg_search(mol, arg) -> function_arg(mol, arg, Search)

#define CTSMARTS_API_ARG_SEARCH(function, arg, Search, search) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig> \
    constexpr auto function##_##arg##_##search(const Molecule::Molecule auto &mol, const auto &arg) \
    { return function##_##arg<SMARTS, Config, SearchType::Search>(mol, arg); }

#define CTSMARTS_API_ATOM_UNIQUE(function) CTSMARTS_API_ARG_SEARCH(function, atom, Unique, unique)
#define CTSMARTS_API_BOND_UNIQUE(function) CTSMARTS_API_ARG_SEARCH(function, bond, Unique, unique)
#define CTSMARTS_API_ATOM_ALL(function) CTSMARTS_API_ARG_SEARCH(function, atom, All, all)
#define CTSMARTS_API_BOND_ALL(function) CTSMARTS_API_ARG_SEARCH(function, bond, All, all)

// caller(mol, arg) -> callee(mol, arg)

#define CTSMARTS_API_OVERLOAD_ARG(caller, callee, Arg, arg) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol> \
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) \
    constexpr auto caller(const Mol &mol, typename Molecule::MoleculeTraits<Mol>::Arg arg) \
    { return callee<SMARTS, Config>(mol, arg); }

#define CTSMARTS_API_OVERLOAD_ATOM(function) CTSMARTS_API_OVERLOAD_ARG(function, function##_atom, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND(function) CTSMARTS_API_OVERLOAD_ARG(function, function##_bond, Bond, bond)

#define CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(function) CTSMARTS_API_OVERLOAD_ARG(function##_unique, function##_atom_##unique, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND_UNIQUE(function) CTSMARTS_API_OVERLOAD_ARG(function##_unique, function##_bond_##unique, Bond, bond)
#define CTSMARTS_API_OVERLOAD_ATOM_ALL(function) CTSMARTS_API_OVERLOAD_ARG(function##_all, function##_atom_##all, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND_ALL(function) CTSMARTS_API_OVERLOAD_ARG(function##_all, function##_bond_##all, Bond, bond)

// function(mol, arg, search) -> function_atom(mol, arg, search)

#define CTSMARTS_API_OVERLOAD_ARG_SEARCH(function, Arg, arg) \
    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol> \
    requires (!Molecule::MoleculeTraits<Mol>::SameAtomBondType) \
    constexpr auto function(const Mol &mol, typename Molecule::MoleculeTraits<Mol>::Arg arg, SearchTypeTag<M> searchType = {}) \
    { return function##_##arg<SMARTS, Config, M>(mol, arg, searchType); }

#define CTSMARTS_API_OVERLOAD_ATOM_SEARCH(function) CTSMARTS_API_OVERLOAD_ARG_SEARCH(function, Atom, atom)
#define CTSMARTS_API_OVERLOAD_BOND_SEARCH(function) CTSMARTS_API_OVERLOAD_ARG_SEARCH(function, Bond, bond)

} // mamespace Kitimar::CTSmarts

 

#include <ctll/list.hpp>

#include <array>

namespace Kitimar::CTSmarts {

    struct Edge
    {
        int index = -1;
        int source = -1;
        int target = -1;

        constexpr auto operator<=>(const Edge&) const noexcept = default;
    };

    namespace impl {

        consteval auto makeEdgeList(auto smarts, auto bonds) noexcept
        {
            if constexpr (ctll::empty(bonds))
                return std::array<Edge, smarts.numBonds>{};
            else {
                auto [bond, tail] = ctll::pop_and_get_front(bonds);
                auto edges = makeEdgeList(smarts, tail);
                edges[bond.index] = Edge{bond.index, bond.source, bond.target};
                return edges;
            }
        }

    } // namespace impl

    template<typename SmartsT>
    struct EdgeList
    {
        static constexpr inline auto data = impl::makeEdgeList(SmartsT{}, SmartsT::bonds);

        consteval EdgeList() noexcept {}
        consteval EdgeList(SmartsT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <array>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto makeVertexDegree(auto smarts, auto edgeList) noexcept
        {
            std::array<int, smarts.numAtoms> degrees = {};
            for (const auto &edge : edgeList.data) {
                ++degrees[edge.source];
                ++degrees[edge.target];
            }
            return degrees;
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT>
    struct VertexDegree
    {
        static constexpr inline auto data = impl::makeVertexDegree(SmartsT{}, EdgeListT{});

        consteval VertexDegree() noexcept {}
        consteval VertexDegree(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <array>
#include <algorithm>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT, typename EdgeListT, typename VertexDegreeT>
        consteval auto makeIncidentList(SmartsT, EdgeListT, VertexDegreeT) noexcept
        {
            constexpr std::size_t stride = *std::ranges::max_element(VertexDegreeT::data);
            std::array<int, SmartsT::numAtoms * stride> incident = {};
            std::ranges::fill(incident, -1);

            std::array<int, SmartsT::numAtoms> sizes = {};
            for (auto i = 0UL; i < SmartsT::numBonds; ++i) {
                auto edge = EdgeListT::data[i];
                auto source = edge.source;
                auto target = edge.target;
                incident[stride * source + sizes[source]] = i;
                incident[stride * target + sizes[target]] = i;
                ++sizes[source];
                ++sizes[target];
            }

            return incident;
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT, typename VertexDegreeT>
    struct IncidentList
    {
        // store adjacent (or incident) bond indices for each vertex
        static constexpr inline auto data = impl::makeIncidentList(SmartsT{}, EdgeListT{}, VertexDegreeT{});
        static constexpr inline auto edges = EdgeListT{};
        static constexpr inline auto degrees = VertexDegreeT{};
        static constexpr inline auto stride = *std::ranges::max_element(VertexDegreeT::data);

        static consteval auto get(int VertexIndex, int IncidentIndex) noexcept
        {
            return data[stride * VertexIndex + IncidentIndex];
        }

        consteval IncidentList() noexcept {}
        consteval IncidentList(SmartsT, EdgeListT, VertexDegreeT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <array>
#include <concepts>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<int SourceIndex, int AdjIndex>
        consteval void dfsSearch(auto smarts, auto incidentList, auto &visitor,
                                 auto &visitedVertices, auto &visitedEdges) noexcept
        {
            constexpr auto sourceDegree = incidentList.degrees.data[SourceIndex];
            if constexpr (AdjIndex < sourceDegree) {
                constexpr auto edgeIndex = incidentList.get(SourceIndex, AdjIndex);

                if (visitedEdges[edgeIndex]) {
                    dfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, visitedVertices, visitedEdges);
                    return;
                }

                constexpr auto edge = incidentList.edges.data[edgeIndex];
                constexpr auto targetIndex = edge.source == SourceIndex ? edge.target : edge.source;
                auto isNewComponent = !visitedVertices[SourceIndex];
                auto isClosure = visitedVertices[targetIndex];

                visitor.visit(edgeIndex, SourceIndex, targetIndex, isNewComponent, isClosure);

                visitedEdges[edgeIndex] = true;
                visitedVertices[SourceIndex] = true;
                visitedVertices[targetIndex] = true;

                // dfs
                if (!isClosure)
                    dfsSearch<targetIndex, 0>(smarts, incidentList, visitor, visitedVertices, visitedEdges);

                visitor.backtrack(edgeIndex, targetIndex, isClosure);

                // next incident bond
                dfsSearch<SourceIndex, AdjIndex + 1>(smarts, incidentList, visitor, visitedVertices, visitedEdges);
            } else if constexpr (SourceIndex == 0) {
                visitor.backtrack(SourceIndex);
            }
        }

    } // namespace impl

    struct DFSVisitorBase
    {
        constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept {}
        constexpr void backtrack(int edge, int target, bool isClosure) noexcept {}
        constexpr void backtrack(int source) noexcept {}
    };

    template<int SourceIndex = 0>
    consteval void dfsSearch(auto smarts, auto incidentList, auto &visitor) noexcept
    {
        std::array<bool, smarts.numAtoms> visitedAtoms = {};
        std::array<bool, smarts.numBonds> visitedBonds = {};
        impl::dfsSearch<SourceIndex, 0>(smarts, incidentList, visitor, visitedAtoms, visitedBonds);
    }

} // namespace Kitimar::CTSmarts

#include <array>

#ifdef KITIMAR_WITH_IOSTREAM
#include <iostream>
#endif

namespace Kitimar::CTSmarts {

    namespace impl {

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

            #ifdef KITIMAR_WITH_IOSTREAM
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
            #endif // KITIMAR_WITH_IOSTREAM

            using Events = std::array<Event, 2 * (SmartsT::numAtoms + SmartsT::numBonds)>;

            consteval DfsSearchEventsVisitor() noexcept {}
            consteval DfsSearchEventsVisitor(SmartsT) noexcept {}

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

        consteval auto makeDfsSearchEvents(auto smarts, auto incidentList) noexcept
        {
            DfsSearchEventsVisitor visitor{smarts};
            dfsSearch(smarts, incidentList, visitor);
            return visitor.events;
        }

    } // namespace impl

    template<typename SmartsT, typename IncidentListT>
    struct DfsSearchEvents
    {
        static constexpr inline auto events = impl::makeDfsSearchEvents(SmartsT{}, IncidentListT{});

        consteval DfsSearchEvents() noexcept {}
        consteval DfsSearchEvents(SmartsT, IncidentListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <ctll/list.hpp>

#include <tuple>
#include <stdexcept>

namespace Kitimar::CTSmarts {

    inline constexpr void PRE(auto cond)
    {
        if (!cond) throw std::runtime_error("");
    }

    template<int Value>
    struct Char : std::integral_constant<char, Value> {};

    // FIXME: rename to Integer
    template<int Value>
    struct Number : std::integral_constant<int, Value> {};

    namespace impl {

        template<typename T>
        consteval T abs(T value) noexcept
        {
            return value >= 0 ? value : -value;
        }

        template<int Coefficient, int Exponent>
        consteval double realValue() noexcept
        {
            auto n = abs(Exponent);
            double mult = 1;
            for (auto i = 0; i < n; ++i)
                mult *= 10;
            return Exponent > 0 ? Coefficient * mult : Coefficient / mult;
        }

    } // namespace impl

    template<int Coefficient, int Exponent = 0>
    struct Real
    {
        static constexpr inline double value = impl::realValue<Coefficient, Exponent>();

        /*
        static consteval double value() noexcept
        {
            auto n = impl::abs(Exponent);
            double mult = 1;
            for (auto i = 0; i < n; ++i)
                mult *= 10;
            return Exponent > 0 ? Coefficient * mult : Coefficient / mult;
        }
        */

        static consteval bool isNear(double n, double tol = 10e-6) noexcept
        {
            //auto delta = value() - n;
            auto delta = value - n;
            return delta * delta < tol * tol;
        }
    };

    //
    // Helper functions to create runtime variable from compile time type in ctll::list
    // See: https://www.scs.stanford.edu/~dm/blog/param-pack.html#array-of-function-pointers
    //

    namespace impl {

        template<std::size_t I, typename R, typename F>
        inline constexpr R with_integral_constant(F f)
        {
            return static_cast<F>(f)(std::integral_constant<std::size_t, I>{});
        }

    } // namespace detail

    template<std::size_t N, typename R = void, typename F>
    inline constexpr R with_n(int n, F &&f)
    {
        constexpr auto invokeArray = [] <std::size_t...I> (std::index_sequence<I...>) {
            return std::array{ impl::with_integral_constant<I, R, F&&>... };
        }(std::make_index_sequence<N>{});

        return invokeArray.at(n)(std::forward<F>(f));
    }

} // namespace ctsmarts

namespace Kitimar::CTSmarts {

    struct DfsEdge : Edge
    {
        bool closure = false;
    };

    namespace impl {

        template<typename SmartsT>
        struct DfsEdgeListVisitor : DFSVisitorBase
        {
            consteval DfsEdgeListVisitor() noexcept {}
            consteval DfsEdgeListVisitor(SmartsT) noexcept {}

            constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
            {
                edges[nextEdgeIndex++] = DfsEdge{{edge, source, target}, isClosure};
            }

            std::array<DfsEdge, SmartsT::numBonds> edges;
            int nextEdgeIndex = 0;
        };

        template<int SourceIndex = 0>
        consteval auto makeDfsEdgeList(auto smarts, auto incidentList) noexcept
        {
            DfsEdgeListVisitor visitor{smarts};
            dfsSearch<SourceIndex>(smarts, incidentList, visitor);
            return visitor.edges;
        }

    } // namespace impl

    template<typename SmartsT, typename IncidentListT, int SourceIndex = 0>
    struct DfsEdgeList
    {
        static constexpr inline auto data = impl::makeDfsEdgeList<SourceIndex>(SmartsT{}, IncidentListT{});

        consteval DfsEdgeList() noexcept {}
        consteval DfsEdgeList(SmartsT, IncidentListT, Number<SourceIndex> = {}) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <array>
#include <tuple>

namespace Kitimar::CTSmarts {

    namespace impl {

        template<typename SmartsT, typename IncidentListT>
        struct CycleMembershipVisitor
        {

            consteval CycleMembershipVisitor() noexcept {}
            consteval CycleMembershipVisitor(SmartsT, IncidentListT) noexcept {}

            constexpr void visit(int edge, int source, int target, bool isNewComponent, bool isClosure) noexcept
            {
                path[depth++] = edge;

                if (isClosure) {
                    for (auto i = depth - 1; i >= 0; --i) {
                        auto e = IncidentListT::edges.data[path[i]];
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

        consteval auto makeCycleMembership(auto smarts, auto incidentList) noexcept
        {
            CycleMembershipVisitor visitor{smarts, incidentList};
            dfsSearch(smarts, incidentList, visitor);
            return std::make_tuple(visitor.vertices, visitor.edges);
        }

    } // namespace impl

    template<typename SmartsT, typename EdgeListT>
    struct CycleMembership
    {
        static constexpr inline auto data = impl::makeCycleMembership(SmartsT{}, EdgeListT{});
        static constexpr inline auto vertices = std::get<0>(data);
        static constexpr inline auto edges = std::get<1>(data);

        consteval CycleMembership() noexcept {}
        consteval CycleMembership(SmartsT, EdgeListT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    template<typename Source, typename Target, typename Expr, bool IsCyclic, bool IsRingClosure>
    struct DfsBond
    {
        static constexpr inline auto source = Source{};
        static constexpr inline auto target = Target{};
        static constexpr inline auto expr = Expr();
        static constexpr inline auto isCyclic = IsCyclic;
        static constexpr inline auto isRingClosure = IsRingClosure;
    };

    namespace impl {

        template<int DfsEdgeIndex>
        consteval auto makeDfsBondListHelper(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            if constexpr (DfsEdgeIndex == smarts.numBonds) {
                return ctll::empty_list{};
            } else {
                constexpr auto edge = dfsEdgeList.data[DfsEdgeIndex];
                constexpr auto expr = get<edge.index>(smarts.bonds).expr;
                auto source = get<edge.source>(smarts.atoms);
                auto target = get<edge.target>(smarts.atoms);
                constexpr auto isCyclic = cycleMembership.edges[edge.index];

                auto dfsBond = DfsBond<decltype(source), decltype(target), decltype(expr), isCyclic, edge.closure>{};

                return ctll::push_front(dfsBond, makeDfsBondListHelper<DfsEdgeIndex + 1>(smarts, dfsEdgeList, cycleMembership));
            }
        }

        consteval auto makeDfsBondList(auto smarts, auto dfsEdgeList, auto cycleMembership) noexcept
        {
            return makeDfsBondListHelper<0>(smarts, dfsEdgeList, cycleMembership);
        }

    } // namespace impl

    template<typename SmartsT, typename DfsEdgeListT, typename CycleMembershipT>
    struct DfsBondList
    {
        static constexpr inline auto data = impl::makeDfsBondList(SmartsT{}, DfsEdgeListT{}, CycleMembershipT{});

        consteval DfsBondList() noexcept {}
        consteval DfsBondList(SmartsT, DfsEdgeListT, CycleMembershipT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Atom primitives
    //

    // '*'
    struct AnyAtom {};

    // 'A'
    struct AnyAliphatic {};

    // 'a'
    struct AnyAromatic {};

    // 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I' | ...
    template<int Element>
    struct AliphaticAtom {};

    // 'b' | 'c' | 'n' | 'o' | 's' | 'p' | 'se' | 'as'
    template<int Element>
    struct AromaticAtom {};

    // NUMBER
    template<int Mass>
    struct Isotope {};

    // '#' NUMBER
    template<int AtomicNumber>
    struct Element {};

    // 'D' | 'D' NUMBER
    template<int N>
    struct Degree {};

    // 'v' | 'v' NUMBER
    template<int N>
    struct Valence {};

    // 'X' | 'X' NUMBER
    template<int N>
    struct Connectivity {};

    // 'H' | 'H' DIGIT
    template<int N>
    struct TotalH {};

    // 'h'
    struct HasImplicitH {};

    // 'h' | 'h' DIGIT
    template<int N>
    struct ImplicitH {};

    // 'R' | 'r' | 'x'
    struct Cyclic {};

    // 'R0' | 'r0' | 'x0'
    struct Acyclic {};

    // 'R' NUMBER
    template<int N>
    struct RingCount {};

    // 'r' NUMBER
    template<int N>
    struct RingSize {};

    // 'x' NUMBER
    template<int N>
    struct RingConnectivity {};

    // '-' | '-' DIGIT | '+' | '+' DIGIT
    template<int N>
    struct Charge
    {
        static constexpr inline int value = N;
    };

    // '@' '?'? | '@@' '?' |
    // '@TH1' '?'? | '@TH2' '?'? |
    // '@AL1' '?'? | '@AL2' '?'? |
    // '@SP1' '?'? | '@SP2' '?'? | '@SP3' '?'? |
    // '@TB1' '?'? | '@TB2' '?'? | '@TB3' '?'? | ... | '@TB20' '?'?
    // '@OH1' '?'? | '@OH2' '?'? | '@OH3' '?'? | ... | '@OH30' '?'?
    // '@TH?' | '@SP?' | '@AL?' | '@TB?' | '@OH?'
    template<int N>
    struct Chiral {};

    // ':' NUMBER
    template<int N>
    struct Class {};

    //
    // Bond primitives
    //

    // single/aromatic depending on context
    struct ImplicitBond {};

    // '-' | '=' | '#' | '$'
    template<int Order>
    struct BondOrder {};

    // '~'
    struct AnyBond {};

    // ':'
    struct AromaticBond {};

    // '@'
    struct RingBond {};

    // '/'
    struct UpBond {};

    // '\'
    struct DownBond {};

    // '/?' | '\?'
    struct UpOrDownBond {};

} // namespace Kitimar::CTSmarts

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    template<typename Expr>
    struct Not
    {
        static constexpr inline auto expr = Expr{};

        constexpr Not() noexcept {}
        constexpr Not(Expr) noexcept {}
    };

    template<typename ...Expr>
    struct And
    {
        static constexpr inline auto expr = ctll::list<Expr...>{};

        constexpr And() noexcept {}
        constexpr And(Expr...) noexcept {}
        constexpr And(ctll::list<Expr...>) noexcept {}
    };

    template<typename ...Expr>
    struct Or
    {
        static constexpr inline auto expr = ctll::list<Expr...>{};

        constexpr Or() noexcept {}
        constexpr Or(Expr...) noexcept {}
        constexpr Or(ctll::list<Expr...>) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    struct AtomExprTag {};
    struct BondExprTag {};

    template<int Index, typename Expr>
    struct Atom
    {
        static constexpr inline auto index = Index;
        static constexpr inline auto expr = Expr{};
    };

    template<int Index, int Source, int Target, typename Expr>
    struct Bond
    {
        static constexpr inline auto index = Index;
        static constexpr inline auto source = Source;
        static constexpr inline auto target = Target;
        static constexpr inline auto expr = Expr{};
    };

    template<typename Atoms, typename Bonds, typename Classes>
    struct BasicSmarts
    {
        static constexpr inline auto atoms = Atoms{};
        static constexpr inline auto bonds = Bonds{};
        static constexpr inline auto classes = Classes{};
        static constexpr inline auto numAtoms = ctll::size(atoms);
        static constexpr inline auto numBonds = ctll::size(bonds);
        static constexpr inline auto numClasses = ctll::size(classes);

        static constexpr inline auto isSingleAtom = numAtoms == 1;
        static constexpr inline auto isSingleBond = numAtoms == 2 && numBonds == 1;

        consteval BasicSmarts() noexcept = default;
        consteval BasicSmarts(Atoms, Bonds, Classes) noexcept {}
    };

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    namespace impl {

        // Primitive
        template<typename Query, typename Expr>
        constexpr auto containsExprHelper(Query, Expr) noexcept
        {
            return std::is_same_v<Query, Expr>;
        }

        template<typename Query, typename ...Expr> constexpr auto containsExprHelper(Query, ctll::list<Expr...>) noexcept;
        template<typename Query, typename ...Expr> constexpr auto containsExprHelper(Query, Or<Expr...>) noexcept;
        template<typename Query, typename ...Expr> constexpr auto containsExprHelper(Query, And<Expr...>) noexcept;

        // Not
        template<typename Query, typename Expr>
        constexpr auto containsExprHelper(Query, Not<Expr>) noexcept
        {
            if constexpr (std::is_same_v<Query, Not<Expr>>)
                return true;
            else
                return containsExprHelper(Query{}, Expr{});
        }

        // Or
        template<typename Query, typename ...Expr>
        constexpr auto containsExprHelper(Query, Or<Expr...>) noexcept
        {
            if constexpr (std::is_same_v<Query, Or<Expr...>>)
                return true;
            else
                return (containsExpr(Query{}, Expr{}) || ...);
        }

        // And
        template<typename Query, typename ...Expr>
        constexpr auto containsExprHelper(Query, And<Expr...>) noexcept
        {
            if constexpr (std::is_same_v<Query, And<Expr...>>)
                return true;
            else
                return (containsExpr(Query{}, Expr{}) || ...);
        }

        // Atom
        template<typename Query, int Index, typename Expr>
        constexpr auto containsExprHelper(Query, Atom<Index, Expr>) noexcept
        {
            return containsExprHelper(Query{}, Expr{});
        }

        // ctll::list
        template<typename Query, typename ...Expr>
        constexpr auto containsExprHelper(Query, ctll::list<Expr...>) noexcept
        {
            return (containsExpr(Query{}, Expr{}) || ...);
        }

    } // namespace impl

    template<typename Query, typename Expr>
    constexpr auto containsExpr(Query, Expr) noexcept
    {
        return impl::containsExprHelper(Query{}, Expr{});
    }

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    namespace impl {

        // Primitive
        template<typename OldExpr, typename NewExpr, typename Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Expr) noexcept
        {
            if constexpr (std::is_same_v<Expr, OldExpr>)
                return NewExpr{};
            else
                return Expr{};
        }

        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, ctll::list<Expr...>) noexcept;
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Or<Expr...>) noexcept;
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, And<Expr...>) noexcept;

        // Not
        template<typename OldExpr, typename NewExpr, typename Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Not<Expr>) noexcept
        {
            if constexpr (std::is_same_v<Not<Expr>, OldExpr>)
                return NewExpr{};
            else
                return Not{replaceExprHelper(OldExpr{}, NewExpr{}, Expr{})};
        }

        // Or
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Or<Expr...>) noexcept
        {
            if constexpr (std::is_same_v<Or<Expr...>, OldExpr>)
                return NewExpr{};
            else
                return Or{replaceExprHelper(OldExpr{}, NewExpr{}, ctll::list<Expr...>{})};
        }

        // And
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, And<Expr...>) noexcept
        {
            if constexpr (std::is_same_v<And<Expr...>, OldExpr>)
                return NewExpr{};
            else
                return And{replaceExprHelper(OldExpr{}, NewExpr{}, ctll::list<Expr...>{})};
        }

        // Atom
        template<typename OldExpr, typename NewExpr, int Index, typename Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, Atom<Index, Expr>) noexcept
        {
            return Atom<Index, decltype(replaceExprHelper(OldExpr{}, NewExpr{}, Expr{}))>{};
        }

        // ctll::list
        template<typename OldExpr, typename NewExpr, typename ...Expr>
        constexpr auto replaceExprHelper(OldExpr, NewExpr, ctll::list<Expr...> l) noexcept
        {
            return decltype(transform(l, [] (auto expr) { return replaceExprHelper(OldExpr{}, NewExpr{}, expr); })){};
        }

    } // namespace impl

    template<typename OldExpr, typename NewExpr, typename Expr>
    constexpr auto replaceExpr(OldExpr, NewExpr, Expr) noexcept
    {
        return impl::replaceExprHelper(OldExpr{}, NewExpr{}, Expr{});
    }

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    // Pubchem (July 2023)
    //
    // # molecules: 111,285,339
    // # atoms:     5,847,611,673
    // # bonds:     6,061,167,923

    // Atom Primitives

    consteval double expressionFrequency(AliphaticAtom<1>) noexcept { return 0.4749471993962647; }
    consteval double expressionFrequency(AliphaticAtom<2>) noexcept { return 2.171826532358993e-08; }
    consteval double expressionFrequency(AliphaticAtom<3>) noexcept { return 1.3528431157263585e-05; }
    consteval double expressionFrequency(AliphaticAtom<4>) noexcept { return 4.076877856309412e-07; }
    consteval double expressionFrequency(AliphaticAtom<5>) noexcept { return 0.00014672845452719994; }
    consteval double expressionFrequency(AliphaticAtom<6>) noexcept { return 0.1961122580168576; }
    consteval double expressionFrequency(AliphaticAtom<7>) noexcept { return 0.03316901603142981; }
    consteval double expressionFrequency(AliphaticAtom<8>) noexcept { return 0.05408450409410525; }
    consteval double expressionFrequency(AliphaticAtom<9>) noexcept { return 0.00826219999340956; }
    consteval double expressionFrequency(AliphaticAtom<10>) noexcept { return 1.0602609452858378e-08; }
    consteval double expressionFrequency(AliphaticAtom<11>) noexcept { return 6.614528128412574e-05; }
    consteval double expressionFrequency(AliphaticAtom<12>) noexcept { return 6.165078608685806e-06; }
    consteval double expressionFrequency(AliphaticAtom<13>) noexcept { return 8.761525270386841e-06; }
    consteval double expressionFrequency(AliphaticAtom<14>) noexcept { return 0.0003513167591350574; }
    consteval double expressionFrequency(AliphaticAtom<15>) noexcept { return 0.0004908912037616586; }
    consteval double expressionFrequency(AliphaticAtom<16>) noexcept { return 0.00477802555367341; }
    consteval double expressionFrequency(AliphaticAtom<17>) noexcept { return 0.004457371105231116; }
    consteval double expressionFrequency(AliphaticAtom<18>) noexcept { return 8.211901995866846e-07; }
    consteval double expressionFrequency(AliphaticAtom<19>) noexcept { return 2.1197550127884277e-05; }
    consteval double expressionFrequency(AliphaticAtom<20>) noexcept { return 4.733213188631183e-06; }
    consteval double expressionFrequency(AliphaticAtom<21>) noexcept { return 3.726307708581941e-07; }
    consteval double expressionFrequency(AliphaticAtom<22>) noexcept { return 8.340499399160852e-06; }
    consteval double expressionFrequency(AliphaticAtom<23>) noexcept { return 6.662547437561399e-06; }
    consteval double expressionFrequency(AliphaticAtom<24>) noexcept { return 3.1927566989965304e-06; }
    consteval double expressionFrequency(AliphaticAtom<25>) noexcept { return 3.548454910407224e-06; }
    consteval double expressionFrequency(AliphaticAtom<26>) noexcept { return 1.1315044241420266e-05; }
    consteval double expressionFrequency(AliphaticAtom<27>) noexcept { return 6.1840621816283525e-06; }
    consteval double expressionFrequency(AliphaticAtom<28>) noexcept { return 6.170726776109029e-06; }
    consteval double expressionFrequency(AliphaticAtom<29>) noexcept { return 9.285843987126784e-06; }
    consteval double expressionFrequency(AliphaticAtom<30>) noexcept { return 9.073104652362688e-06; }
    consteval double expressionFrequency(AliphaticAtom<31>) noexcept { return 1.402110691959276e-06; }
    consteval double expressionFrequency(AliphaticAtom<32>) noexcept { return 6.8320210378386795e-06; }
    consteval double expressionFrequency(AliphaticAtom<33>) noexcept { return 4.720216693008036e-06; }
    consteval double expressionFrequency(AliphaticAtom<34>) noexcept { return 1.450626448327788e-05; }
    consteval double expressionFrequency(AliphaticAtom<35>) noexcept { return 0.0015151734741764703; }
    consteval double expressionFrequency(AliphaticAtom<36>) noexcept { return 1.5561915997712353e-08; }
    consteval double expressionFrequency(AliphaticAtom<37>) noexcept { return 2.469555387640306e-06; }
    consteval double expressionFrequency(AliphaticAtom<38>) noexcept { return 6.717270759248571e-07; }
    consteval double expressionFrequency(AliphaticAtom<39>) noexcept { return 2.6826331493392595e-05; }
    consteval double expressionFrequency(AliphaticAtom<40>) noexcept { return 1.1579084791229487e-05; }
    consteval double expressionFrequency(AliphaticAtom<41>) noexcept { return 5.383393050033839e-07; }
    consteval double expressionFrequency(AliphaticAtom<42>) noexcept { return 2.751036850827957e-06; }
    consteval double expressionFrequency(AliphaticAtom<43>) noexcept { return 5.778426455255303e-07; }
    consteval double expressionFrequency(AliphaticAtom<44>) noexcept { return 5.253596889215508e-06; }
    consteval double expressionFrequency(AliphaticAtom<45>) noexcept { return 2.2660526999497563e-06; }
    consteval double expressionFrequency(AliphaticAtom<46>) noexcept { return 6.156017290436009e-06; }
    consteval double expressionFrequency(AliphaticAtom<47>) noexcept { return 2.9815587809623268e-06; }
    consteval double expressionFrequency(AliphaticAtom<48>) noexcept { return 1.0896754338453948e-06; }
    consteval double expressionFrequency(AliphaticAtom<49>) noexcept { return 1.2815487308949387e-06; }
    consteval double expressionFrequency(AliphaticAtom<50>) noexcept { return 1.5540705027416677e-05; }
    consteval double expressionFrequency(AliphaticAtom<51>) noexcept { return 2.518804315376452e-06; }
    consteval double expressionFrequency(AliphaticAtom<52>) noexcept { return 3.027217542547342e-06; }
    consteval double expressionFrequency(AliphaticAtom<53>) noexcept { return 0.00044532839948324926; }
    consteval double expressionFrequency(AliphaticAtom<54>) noexcept { return 7.011410800521098e-08; }
    consteval double expressionFrequency(AliphaticAtom<55>) noexcept { return 1.2071586531912752e-06; }
    consteval double expressionFrequency(AliphaticAtom<56>) noexcept { return 1.6416964225963227e-06; }
    consteval double expressionFrequency(AliphaticAtom<57>) noexcept { return 6.08453612724821e-07; }
    consteval double expressionFrequency(AliphaticAtom<58>) noexcept { return 6.821588400097205e-07; }
    consteval double expressionFrequency(AliphaticAtom<59>) noexcept { return 4.748948343291471e-07; }
    consteval double expressionFrequency(AliphaticAtom<60>) noexcept { return 5.143983651460342e-07; }
    consteval double expressionFrequency(AliphaticAtom<61>) noexcept { return 3.676714616283999e-08; }
    consteval double expressionFrequency(AliphaticAtom<62>) noexcept { return 4.365885665702247e-07; }
    consteval double expressionFrequency(AliphaticAtom<63>) noexcept { return 6.180301544863463e-07; }
    consteval double expressionFrequency(AliphaticAtom<64>) noexcept { return 1.0932663098649765e-06; }
    consteval double expressionFrequency(AliphaticAtom<65>) noexcept { return 1.7930371706930948e-06; }
    consteval double expressionFrequency(AliphaticAtom<66>) noexcept { return 3.669873997759469e-07; }
    consteval double expressionFrequency(AliphaticAtom<67>) noexcept { return 2.0965824610991902e-07; }
    consteval double expressionFrequency(AliphaticAtom<68>) noexcept { return 2.633553715007652e-07; }
    consteval double expressionFrequency(AliphaticAtom<69>) noexcept { return 1.4570050273142341e-07; }
    consteval double expressionFrequency(AliphaticAtom<70>) noexcept { return 4.1076589223397823e-07; }
    consteval double expressionFrequency(AliphaticAtom<71>) noexcept { return 2.6797268759377586e-07; }
    consteval double expressionFrequency(AliphaticAtom<72>) noexcept { return 2.031769833345076e-06; }
    consteval double expressionFrequency(AliphaticAtom<73>) noexcept { return 4.899436548678864e-07; }
    consteval double expressionFrequency(AliphaticAtom<74>) noexcept { return 1.001365986837606e-05; }
    consteval double expressionFrequency(AliphaticAtom<75>) noexcept { return 1.441784550958852e-06; }
    consteval double expressionFrequency(AliphaticAtom<76>) noexcept { return 1.3672254091369582e-06; }
    consteval double expressionFrequency(AliphaticAtom<77>) noexcept { return 1.9995865327071727e-05; }
    consteval double expressionFrequency(AliphaticAtom<78>) noexcept { return 1.7766228466835344e-05; }
    consteval double expressionFrequency(AliphaticAtom<79>) noexcept { return 2.7898568890259374e-06; }
    consteval double expressionFrequency(AliphaticAtom<80>) noexcept { return 2.427143445630615e-06; }
    consteval double expressionFrequency(AliphaticAtom<81>) noexcept { return 9.268742575011324e-07; }
    consteval double expressionFrequency(AliphaticAtom<82>) noexcept { return 1.9088139148525278e-06; }
    consteval double expressionFrequency(AliphaticAtom<83>) noexcept { return 1.5088208143517404e-06; }
    consteval double expressionFrequency(AliphaticAtom<84>) noexcept { return 1.6656372995672705e-07; }
    consteval double expressionFrequency(AliphaticAtom<85>) noexcept { return 9.67916525511301e-08; }
    consteval double expressionFrequency(AliphaticAtom<86>) noexcept { return 6.669395624367344e-09; }
    consteval double expressionFrequency(AliphaticAtom<87>) noexcept { return 2.9071689765771146e-09; }
    consteval double expressionFrequency(AliphaticAtom<88>) noexcept { return 3.557007235641863e-08; }
    consteval double expressionFrequency(AliphaticAtom<89>) noexcept { return 2.2439941712741106e-06; }
    consteval double expressionFrequency(AliphaticAtom<90>) noexcept { return 2.0025271690284677e-07; }
    consteval double expressionFrequency(AliphaticAtom<91>) noexcept { return 2.411241145565803e-08; }
    consteval double expressionFrequency(AliphaticAtom<92>) noexcept { return 3.4993779409602786e-06; }
    consteval double expressionFrequency(AliphaticAtom<93>) noexcept { return 6.960108567889844e-08; }
    consteval double expressionFrequency(AliphaticAtom<94>) noexcept { return 5.2329022818282504e-08; }
    consteval double expressionFrequency(AliphaticAtom<95>) noexcept { return 4.15554213571829e-08; }
    consteval double expressionFrequency(AliphaticAtom<96>) noexcept { return 7.028510775364663e-08; }
    consteval double expressionFrequency(AliphaticAtom<97>) noexcept { return 1.4022822112723246e-08; }
    consteval double expressionFrequency(AliphaticAtom<98>) noexcept { return 3.1123839878315475e-08; }
    consteval double expressionFrequency(AliphaticAtom<99>) noexcept { return 3.385994289574164e-08; }
    consteval double expressionFrequency(AliphaticAtom<100>) noexcept { return 2.0880316907909633e-07; }
    consteval double expressionFrequency(AliphaticAtom<101>) noexcept { return 6.498373746262715e-09; }
    consteval double expressionFrequency(AliphaticAtom<102>) noexcept { return 7.011409771494334e-08; }
    consteval double expressionFrequency(AliphaticAtom<103>) noexcept { return 6.703591137699363e-08; }
    consteval double expressionFrequency(AliphaticAtom<104>) noexcept { return 1.3472174004981822e-06; }
    consteval double expressionFrequency(Element<1>) noexcept { return 0.4749471993962647; }
    consteval double expressionFrequency(Element<2>) noexcept { return 2.171826532358993e-08; }
    consteval double expressionFrequency(Element<3>) noexcept { return 1.3528431157263585e-05; }
    consteval double expressionFrequency(Element<4>) noexcept { return 4.076877856309412e-07; }
    consteval double expressionFrequency(Element<5>) noexcept { return 0.00017143225264344694; }
    consteval double expressionFrequency(Element<6>) noexcept { return 0.39351618247954073; }
    consteval double expressionFrequency(Element<7>) noexcept { return 0.05256898187077755; }
    consteval double expressionFrequency(Element<8>) noexcept { return 0.05580677409443103; }
    consteval double expressionFrequency(Element<9>) noexcept { return 0.00826219999340956; }
    consteval double expressionFrequency(Element<10>) noexcept { return 1.0602609452858378e-08; }
    consteval double expressionFrequency(Element<11>) noexcept { return 6.614528128412574e-05; }
    consteval double expressionFrequency(Element<12>) noexcept { return 6.165078608685806e-06; }
    consteval double expressionFrequency(Element<13>) noexcept { return 8.761525270386841e-06; }
    consteval double expressionFrequency(Element<14>) noexcept { return 0.0003513167591350574; }
    consteval double expressionFrequency(Element<15>) noexcept { return 0.0004951671528170457; }
    consteval double expressionFrequency(Element<16>) noexcept { return 0.00707726374200449; }
    consteval double expressionFrequency(Element<17>) noexcept { return 0.004457371105231116; }
    consteval double expressionFrequency(Element<18>) noexcept { return 8.211901995866846e-07; }
    consteval double expressionFrequency(Element<19>) noexcept { return 2.1197550127884277e-05; }
    consteval double expressionFrequency(Element<20>) noexcept { return 4.733213188631183e-06; }
    consteval double expressionFrequency(Element<21>) noexcept { return 3.726307708581941e-07; }
    consteval double expressionFrequency(Element<22>) noexcept { return 8.340499399160852e-06; }
    consteval double expressionFrequency(Element<23>) noexcept { return 6.662547437561399e-06; }
    consteval double expressionFrequency(Element<24>) noexcept { return 3.1927566989965304e-06; }
    consteval double expressionFrequency(Element<25>) noexcept { return 3.548454910407224e-06; }
    consteval double expressionFrequency(Element<26>) noexcept { return 1.1315044241420266e-05; }
    consteval double expressionFrequency(Element<27>) noexcept { return 6.1840621816283525e-06; }
    consteval double expressionFrequency(Element<28>) noexcept { return 6.170726776109029e-06; }
    consteval double expressionFrequency(Element<29>) noexcept { return 9.285843987126784e-06; }
    consteval double expressionFrequency(Element<30>) noexcept { return 9.073104652362688e-06; }
    consteval double expressionFrequency(Element<31>) noexcept { return 1.402110691959276e-06; }
    consteval double expressionFrequency(Element<32>) noexcept { return 6.8320210378386795e-06; }
    consteval double expressionFrequency(Element<33>) noexcept { return 4.822139527103718e-06; }
    consteval double expressionFrequency(Element<34>) noexcept { return 1.9262057011545827e-05; }
    consteval double expressionFrequency(Element<35>) noexcept { return 0.0015151734741764703; }
    consteval double expressionFrequency(Element<36>) noexcept { return 1.5561915997712353e-08; }
    consteval double expressionFrequency(Element<37>) noexcept { return 2.469555387640306e-06; }
    consteval double expressionFrequency(Element<38>) noexcept { return 6.717270759248571e-07; }
    consteval double expressionFrequency(Element<39>) noexcept { return 2.6826331493392595e-05; }
    consteval double expressionFrequency(Element<40>) noexcept { return 1.1579084791229487e-05; }
    consteval double expressionFrequency(Element<41>) noexcept { return 5.383393050033839e-07; }
    consteval double expressionFrequency(Element<42>) noexcept { return 2.751036850827957e-06; }
    consteval double expressionFrequency(Element<43>) noexcept { return 5.778426455255303e-07; }
    consteval double expressionFrequency(Element<44>) noexcept { return 5.253596889215508e-06; }
    consteval double expressionFrequency(Element<45>) noexcept { return 2.2660526999497563e-06; }
    consteval double expressionFrequency(Element<46>) noexcept { return 6.156017290436009e-06; }
    consteval double expressionFrequency(Element<47>) noexcept { return 2.9815587809623268e-06; }
    consteval double expressionFrequency(Element<48>) noexcept { return 1.0896754338453948e-06; }
    consteval double expressionFrequency(Element<49>) noexcept { return 1.2815487308949387e-06; }
    consteval double expressionFrequency(Element<50>) noexcept { return 1.5540705027416677e-05; }
    consteval double expressionFrequency(Element<51>) noexcept { return 2.518804315376452e-06; }
    consteval double expressionFrequency(Element<52>) noexcept { return 3.027217542547342e-06; }
    consteval double expressionFrequency(Element<53>) noexcept { return 0.00044532839948324926; }
    consteval double expressionFrequency(Element<54>) noexcept { return 7.011410800521098e-08; }
    consteval double expressionFrequency(Element<55>) noexcept { return 1.2071586531912752e-06; }
    consteval double expressionFrequency(Element<56>) noexcept { return 1.6416964225963227e-06; }
    consteval double expressionFrequency(Element<57>) noexcept { return 6.08453612724821e-07; }
    consteval double expressionFrequency(Element<58>) noexcept { return 6.821588400097205e-07; }
    consteval double expressionFrequency(Element<59>) noexcept { return 4.748948343291471e-07; }
    consteval double expressionFrequency(Element<60>) noexcept { return 5.143983651460342e-07; }
    consteval double expressionFrequency(Element<61>) noexcept { return 3.676714616283999e-08; }
    consteval double expressionFrequency(Element<62>) noexcept { return 4.365885665702247e-07; }
    consteval double expressionFrequency(Element<63>) noexcept { return 6.180301544863463e-07; }
    consteval double expressionFrequency(Element<64>) noexcept { return 1.0932663098649765e-06; }
    consteval double expressionFrequency(Element<65>) noexcept { return 1.7930371706930948e-06; }
    consteval double expressionFrequency(Element<66>) noexcept { return 3.669873997759469e-07; }
    consteval double expressionFrequency(Element<67>) noexcept { return 2.0965824610991902e-07; }
    consteval double expressionFrequency(Element<68>) noexcept { return 2.633553715007652e-07; }
    consteval double expressionFrequency(Element<69>) noexcept { return 1.4570050273142341e-07; }
    consteval double expressionFrequency(Element<70>) noexcept { return 4.1076589223397823e-07; }
    consteval double expressionFrequency(Element<71>) noexcept { return 2.6797268759377586e-07; }
    consteval double expressionFrequency(Element<72>) noexcept { return 2.031769833345076e-06; }
    consteval double expressionFrequency(Element<73>) noexcept { return 4.899436548678864e-07; }
    consteval double expressionFrequency(Element<74>) noexcept { return 1.001365986837606e-05; }
    consteval double expressionFrequency(Element<75>) noexcept { return 1.441784550958852e-06; }
    consteval double expressionFrequency(Element<76>) noexcept { return 1.3672254091369582e-06; }
    consteval double expressionFrequency(Element<77>) noexcept { return 1.9995865327071727e-05; }
    consteval double expressionFrequency(Element<78>) noexcept { return 1.7766228466835344e-05; }
    consteval double expressionFrequency(Element<79>) noexcept { return 2.7898568890259374e-06; }
    consteval double expressionFrequency(Element<80>) noexcept { return 2.427143445630615e-06; }
    consteval double expressionFrequency(Element<81>) noexcept { return 9.268742575011324e-07; }
    consteval double expressionFrequency(Element<82>) noexcept { return 1.9088139148525278e-06; }
    consteval double expressionFrequency(Element<83>) noexcept { return 1.5088208143517404e-06; }
    consteval double expressionFrequency(Element<84>) noexcept { return 1.6656372995672705e-07; }
    consteval double expressionFrequency(Element<85>) noexcept { return 9.67916525511301e-08; }
    consteval double expressionFrequency(Element<86>) noexcept { return 6.669395624367344e-09; }
    consteval double expressionFrequency(Element<87>) noexcept { return 2.9071689765771146e-09; }
    consteval double expressionFrequency(Element<88>) noexcept { return 3.557007235641863e-08; }
    consteval double expressionFrequency(Element<89>) noexcept { return 2.2439941712741106e-06; }
    consteval double expressionFrequency(Element<90>) noexcept { return 2.0025271690284677e-07; }
    consteval double expressionFrequency(Element<91>) noexcept { return 2.411241145565803e-08; }
    consteval double expressionFrequency(Element<92>) noexcept { return 3.4993779409602786e-06; }
    consteval double expressionFrequency(Element<93>) noexcept { return 6.960108567889844e-08; }
    consteval double expressionFrequency(Element<94>) noexcept { return 5.2329022818282504e-08; }
    consteval double expressionFrequency(Element<95>) noexcept { return 4.15554213571829e-08; }
    consteval double expressionFrequency(Element<96>) noexcept { return 7.028510775364663e-08; }
    consteval double expressionFrequency(Element<97>) noexcept { return 1.4022822112723246e-08; }
    consteval double expressionFrequency(Element<98>) noexcept { return 3.1123839878315475e-08; }
    consteval double expressionFrequency(Element<99>) noexcept { return 3.385994289574164e-08; }
    consteval double expressionFrequency(Element<100>) noexcept { return 2.0880316907909633e-07; }
    consteval double expressionFrequency(Element<101>) noexcept { return 6.498373746262715e-09; }
    consteval double expressionFrequency(Element<102>) noexcept { return 7.011409771494334e-08; }
    consteval double expressionFrequency(Element<103>) noexcept { return 6.703591137699363e-08; }
    consteval double expressionFrequency(Element<104>) noexcept { return 1.3472174004981822e-06; }
    consteval double expressionFrequency(Isotope<0>) noexcept { return 0.9996383414418606; }
    consteval double expressionFrequency(Isotope<1>) noexcept { return 7.4218304694171e-08; }
    consteval double expressionFrequency(Isotope<2>) noexcept { return 0.00033513643736927533; }
    consteval double expressionFrequency(Isotope<3>) noexcept { return 8.518005996726257e-06; }
    consteval double expressionFrequency(Isotope<4>) noexcept { return 5.130296707458365e-10; }
    consteval double expressionFrequency(Isotope<5>) noexcept { return 0.0; }
    consteval double expressionFrequency(Isotope<6>) noexcept { return 4.959286874890076e-09; }
    consteval double expressionFrequency(Isotope<7>) noexcept { return 4.617275268425062e-09; }
    consteval double expressionFrequency(Isotope<8>) noexcept { return 1.19707038787466e-09; }
    consteval double expressionFrequency(Isotope<9>) noexcept { return 1.36807994418259e-09; }
    consteval double expressionFrequency(Isotope<10>) noexcept { return 2.204318537121106e-07; }
    consteval double expressionFrequency(Isotope<11>) noexcept { return 6.645445594665464e-07; }
    consteval double expressionFrequency(Isotope<12>) noexcept { return 9.593657226479808e-08; }
    consteval double expressionFrequency(Isotope<13>) noexcept { return 9.56612862675585e-06; }
    consteval double expressionFrequency(Isotope<14>) noexcept { return 1.678805169274481e-06; }
    consteval double expressionFrequency(Isotope<15>) noexcept { return 1.6671763925937938e-06; }
    consteval double expressionFrequency(Isotope<16>) noexcept { return 7.695452647724006e-08; }
    consteval double expressionFrequency(Isotope<17>) noexcept { return 6.464177020227686e-08; }
    consteval double expressionFrequency(Isotope<18>) noexcept { return 2.13676987617894e-06; }
    consteval double expressionFrequency(Isotope<19>) noexcept { return 8.003268942175714e-08; }
    consteval double expressionFrequency(Isotope<20>) noexcept { return 1.368080210426227e-09; }
    consteval double expressionFrequency(Isotope<21>) noexcept { return 3.420200784585e-10; }
    consteval double expressionFrequency(Isotope<22>) noexcept { return 1.5390906375225888e-09; }
    consteval double expressionFrequency(Isotope<23>) noexcept { return 3.420201093524973e-10; }
    consteval double expressionFrequency(Isotope<24>) noexcept { return 1.0260598672207349e-09; }
    consteval double expressionFrequency(Isotope<25>) noexcept { return 2.3941399356015e-09; }
    consteval double expressionFrequency(Isotope<26>) noexcept { return 3.4202027536788187e-10; }
    consteval double expressionFrequency(Isotope<27>) noexcept { return 5.130302671474606e-10; }
    consteval double expressionFrequency(Isotope<28>) noexcept { return 1.710098349194878e-09; }
    consteval double expressionFrequency(Isotope<29>) noexcept { return 5.130302315841724e-10; }
    consteval double expressionFrequency(Isotope<30>) noexcept { return 5.130298350821593e-10; }
    consteval double expressionFrequency(Isotope<31>) noexcept { return 1.8811102968650393e-09; }
    consteval double expressionFrequency(Isotope<32>) noexcept { return 3.214985642204891e-08; }
    consteval double expressionFrequency(Isotope<33>) noexcept { return 5.814346583481931e-09; }
    consteval double expressionFrequency(Isotope<34>) noexcept { return 8.550506714540189e-09; }
    consteval double expressionFrequency(Isotope<35>) noexcept { return 4.343653775524278e-08; }
    consteval double expressionFrequency(Isotope<36>) noexcept { return 3.76221927074334e-09; }
    consteval double expressionFrequency(Isotope<37>) noexcept { return 1.0431607106519997e-08; }
    consteval double expressionFrequency(Isotope<38>) noexcept { return 1.7101009924878782e-09; }
    consteval double expressionFrequency(Isotope<39>) noexcept { return 6.840401272813281e-10; }
    consteval double expressionFrequency(Isotope<40>) noexcept { return 1.1970700287859557e-09; }
    consteval double expressionFrequency(Isotope<41>) noexcept { return 6.840402873034282e-10; }
    consteval double expressionFrequency(Isotope<42>) noexcept { return 8.550493666910602e-10; }
    consteval double expressionFrequency(Isotope<43>) noexcept { return 1.197069757469786e-09; }
    consteval double expressionFrequency(Isotope<44>) noexcept { return 1.539089785771775e-09; }
    consteval double expressionFrequency(Isotope<45>) noexcept { return 1.3680797907947367e-09; }
    consteval double expressionFrequency(Isotope<46>) noexcept { return 5.130299872780215e-10; }
    consteval double expressionFrequency(Isotope<47>) noexcept { return 1.1970704351820253e-09; }
    consteval double expressionFrequency(Isotope<48>) noexcept { return 8.550499428914456e-10; }
    consteval double expressionFrequency(Isotope<49>) noexcept { return 8.550499428914456e-10; }
    consteval double expressionFrequency(Isotope<50>) noexcept { return 5.130302813216065e-10; }
    consteval double expressionFrequency(Isotope<51>) noexcept { return 2.736159652084903e-09; }
    consteval double expressionFrequency(Isotope<52>) noexcept { return 1.3680803462006155e-09; }
    consteval double expressionFrequency(Isotope<53>) noexcept { return 6.840398482094623e-10; }
    consteval double expressionFrequency(Isotope<54>) noexcept { return 5.130299181942789e-10; }
    consteval double expressionFrequency(Isotope<55>) noexcept { return 1.1970698174463918e-09; }
    consteval double expressionFrequency(Isotope<56>) noexcept { return 6.840399819093118e-10; }
    consteval double expressionFrequency(Isotope<57>) noexcept { return 4.104242632690942e-09; }
    consteval double expressionFrequency(Isotope<58>) noexcept { return 1.881110913798157e-09; }
    consteval double expressionFrequency(Isotope<59>) noexcept { return 3.24919086847905e-09; }
    consteval double expressionFrequency(Isotope<60>) noexcept { return 3.762218244128093e-09; }
    consteval double expressionFrequency(Isotope<61>) noexcept { return 8.550497104233412e-10; }
    consteval double expressionFrequency(Isotope<62>) noexcept { return 2.7361597765495977e-09; }
    consteval double expressionFrequency(Isotope<63>) noexcept { return 8.55050231013023e-10; }
    consteval double expressionFrequency(Isotope<64>) noexcept { return 6.241863921878742e-08; }
    consteval double expressionFrequency(Isotope<65>) noexcept { return 1.7100992767986084e-09; }
    consteval double expressionFrequency(Isotope<66>) noexcept { return 1.8811066060065647e-09; }
    consteval double expressionFrequency(Isotope<67>) noexcept { return 1.7956049118596935e-08; }
    consteval double expressionFrequency(Isotope<68>) noexcept { return 3.779319714051467e-08; }
    consteval double expressionFrequency(Isotope<69>) noexcept { return 1.026060707389726e-09; }
    consteval double expressionFrequency(Isotope<70>) noexcept { return 1.7100988068362625e-09; }
    consteval double expressionFrequency(Isotope<71>) noexcept { return 1.0260596542452211e-09; }
    consteval double expressionFrequency(Isotope<72>) noexcept { return 1.5390893097264157e-09; }
    consteval double expressionFrequency(Isotope<73>) noexcept { return 1.0260591346655587e-09; }
    consteval double expressionFrequency(Isotope<74>) noexcept { return 2.565149805019243e-09; }
    consteval double expressionFrequency(Isotope<75>) noexcept { return 8.721505406150352e-09; }
    consteval double expressionFrequency(Isotope<76>) noexcept { return 3.3688960217641203e-08; }
    consteval double expressionFrequency(Isotope<77>) noexcept { return 1.4193818383247572e-08; }
    consteval double expressionFrequency(Isotope<78>) noexcept { return 8.550503173372401e-10; }
    consteval double expressionFrequency(Isotope<79>) noexcept { return 3.0781782482766106e-09; }
    consteval double expressionFrequency(Isotope<80>) noexcept { return 5.9853438409892e-09; }
    consteval double expressionFrequency(Isotope<81>) noexcept { return 2.736161195553811e-09; }
    consteval double expressionFrequency(Isotope<82>) noexcept { return 5.6433313469494475e-09; }
    consteval double expressionFrequency(Isotope<83>) noexcept { return 1.1970702675715793e-09; }
    consteval double expressionFrequency(Isotope<84>) noexcept { return 6.840400575224039e-10; }
    consteval double expressionFrequency(Isotope<85>) noexcept { return 1.7100997635486628e-09; }
    consteval double expressionFrequency(Isotope<86>) noexcept { return 2.0521209024015163e-09; }
    consteval double expressionFrequency(Isotope<87>) noexcept { return 1.5390895221051728e-09; }
    consteval double expressionFrequency(Isotope<88>) noexcept { return 1.368079729475959e-09; }
    consteval double expressionFrequency(Isotope<89>) noexcept { return 5.985349766512155e-09; }
    consteval double expressionFrequency(Isotope<90>) noexcept { return 1.4535851276945658e-08; }
    consteval double expressionFrequency(Isotope<91>) noexcept { return 8.5505015394703e-10; }
    consteval double expressionFrequency(Isotope<92>) noexcept { return 1.0260596939027795e-09; }
    consteval double expressionFrequency(Isotope<93>) noexcept { return 8.550499945325456e-10; }
    consteval double expressionFrequency(Isotope<94>) noexcept { return 1.710099221824746e-09; }
    consteval double expressionFrequency(Isotope<95>) noexcept { return 1.5390900303514563e-09; }
    consteval double expressionFrequency(Isotope<96>) noexcept { return 6.840401800787978e-10; }
    consteval double expressionFrequency(Isotope<97>) noexcept { return 8.550499513489667e-10; }
    consteval double expressionFrequency(Isotope<98>) noexcept { return 5.472318215007247e-09; }
    consteval double expressionFrequency(Isotope<99>) noexcept { return 2.9892557837097385e-07; }
    consteval double expressionFrequency(Isotope<100>) noexcept { return 1.3680796761414495e-09; }
    consteval double expressionFrequency(Isotope<101>) noexcept { return 8.550500622020357e-10; }
    consteval double expressionFrequency(Isotope<102>) noexcept { return 6.840399819093118e-10; }
    consteval double expressionFrequency(Isotope<103>) noexcept { return 1.5390903782367495e-09; }
    consteval double expressionFrequency(Isotope<104>) noexcept { return 8.550498705317329e-10; }
    consteval double expressionFrequency(Isotope<105>) noexcept { return 1.197069057319395e-09; }
    consteval double expressionFrequency(Isotope<106>) noexcept { return 1.3680801345831605e-09; }
    consteval double expressionFrequency(Isotope<107>) noexcept { return 6.840398097340962e-10; }
    consteval double expressionFrequency(Isotope<108>) noexcept { return 5.13030139956679e-10; }
    consteval double expressionFrequency(Isotope<109>) noexcept { return 1.1970699741231603e-09; }
    consteval double expressionFrequency(Isotope<110>) noexcept { return 1.3680796832750456e-09; }
    consteval double expressionFrequency(Isotope<111>) noexcept { return 3.7451177655750345e-08; }
    consteval double expressionFrequency(Isotope<112>) noexcept { return 8.550497539797732e-10; }
    consteval double expressionFrequency(Isotope<113>) noexcept { return 4.7882800061306495e-09; }
    consteval double expressionFrequency(Isotope<114>) noexcept { return 8.550496647954823e-10; }
    consteval double expressionFrequency(Isotope<115>) noexcept { return 1.0260596220303256e-09; }
    consteval double expressionFrequency(Isotope<116>) noexcept { return 8.550509920267409e-10; }
    consteval double expressionFrequency(Isotope<117>) noexcept { return 3.933228686484395e-09; }
    consteval double expressionFrequency(Isotope<118>) noexcept { return 6.840397865101532e-10; }
    consteval double expressionFrequency(Isotope<119>) noexcept { return 1.3680798691871343e-09; }
    consteval double expressionFrequency(Isotope<120>) noexcept { return 1.368079478682014e-09; }
    consteval double expressionFrequency(Isotope<121>) noexcept { return 3.9332292621484235e-09; }
    consteval double expressionFrequency(Isotope<122>) noexcept { return 2.22312926803823e-09; }
    consteval double expressionFrequency(Isotope<123>) noexcept { return 2.1393347139327893e-07; }
    consteval double expressionFrequency(Isotope<124>) noexcept { return 6.139260389805841e-08; }
    consteval double expressionFrequency(Isotope<125>) noexcept { return 3.0593687596924825e-07; }
    consteval double expressionFrequency(Isotope<126>) noexcept { return 1.5390900284116728e-09; }
    consteval double expressionFrequency(Isotope<127>) noexcept { return 6.1563543757339895e-09; }
    consteval double expressionFrequency(Isotope<128>) noexcept { return 1.3680809316680326e-09; }
    consteval double expressionFrequency(Isotope<129>) noexcept { return 2.56514950933887e-09; }
    consteval double expressionFrequency(Isotope<130>) noexcept { return 1.3680801388718686e-09; }
    consteval double expressionFrequency(Isotope<131>) noexcept { return 1.16286780274583e-07; }
    consteval double expressionFrequency(Isotope<132>) noexcept { return 1.5390892301600173e-09; }
    consteval double expressionFrequency(Isotope<133>) noexcept { return 1.1970699861589288e-09; }
    consteval double expressionFrequency(Isotope<134>) noexcept { return 1.5390894417038356e-09; }
    consteval double expressionFrequency(Isotope<135>) noexcept { return 1.5390893282816833e-09; }
    consteval double expressionFrequency(Isotope<136>) noexcept { return 1.3680796254324907e-09; }
    consteval double expressionFrequency(Isotope<137>) noexcept { return 1.7100997723876049e-09; }
    consteval double expressionFrequency(Isotope<138>) noexcept { return 1.3680798529558923e-09; }
    consteval double expressionFrequency(Isotope<139>) noexcept { return 1.0260599328785149e-09; }
    consteval double expressionFrequency(Isotope<140>) noexcept { return 5.130301089534575e-10; }
    consteval double expressionFrequency(Isotope<141>) noexcept { return 1.710099774625527e-09; }
    consteval double expressionFrequency(Isotope<142>) noexcept { return 1.1970703060962818e-09; }
    consteval double expressionFrequency(Isotope<143>) noexcept { return 1.0260600423735934e-09; }
    consteval double expressionFrequency(Isotope<144>) noexcept { return 1.881108968696099e-09; }
    consteval double expressionFrequency(Isotope<145>) noexcept { return 1.026059929155679e-09; }
    consteval double expressionFrequency(Isotope<146>) noexcept { return 1.0260597764770211e-09; }
    consteval double expressionFrequency(Isotope<147>) noexcept { return 1.5390901648280332e-09; }
    consteval double expressionFrequency(Isotope<148>) noexcept { return 1.1970701272184152e-09; }
    consteval double expressionFrequency(Isotope<149>) noexcept { return 2.3941414499827374e-09; }
    consteval double expressionFrequency(Isotope<150>) noexcept { return 8.550500955701031e-10; }
    consteval double expressionFrequency(Isotope<151>) noexcept { return 1.197069834379134e-09; }
    consteval double expressionFrequency(Isotope<152>) noexcept { return 1.5390896330718606e-09; }
    consteval double expressionFrequency(Isotope<153>) noexcept { return 6.6693876902053465e-09; }
    consteval double expressionFrequency(Isotope<154>) noexcept { return 8.550497553337107e-10; }
    consteval double expressionFrequency(Isotope<155>) noexcept { return 1.0260599375227357e-09; }
    consteval double expressionFrequency(Isotope<156>) noexcept { return 1.0260599744934876e-09; }
    consteval double expressionFrequency(Isotope<157>) noexcept { return 1.3680809924362433e-09; }
    consteval double expressionFrequency(Isotope<158>) noexcept { return 6.840401177633911e-10; }
    consteval double expressionFrequency(Isotope<159>) noexcept { return 1.1970698886432708e-09; }
    consteval double expressionFrequency(Isotope<160>) noexcept { return 5.13030139956679e-10; }
    consteval double expressionFrequency(Isotope<161>) noexcept { return 1.0260596332113005e-09; }
    consteval double expressionFrequency(Isotope<162>) noexcept { return 8.550500955701031e-10; }
    consteval double expressionFrequency(Isotope<163>) noexcept { return 1.7101018434325505e-10; }
    consteval double expressionFrequency(Isotope<164>) noexcept { return 1.0260600837881768e-09; }
    consteval double expressionFrequency(Isotope<165>) noexcept { return 1.3680793215124974e-09; }
    consteval double expressionFrequency(Isotope<166>) noexcept { return 4.275248467659388e-09; }
    consteval double expressionFrequency(Isotope<167>) noexcept { return 1.0260599471735646e-09; }
    consteval double expressionFrequency(Isotope<168>) noexcept { return 5.130297997202866e-10; }
    consteval double expressionFrequency(Isotope<169>) noexcept { return 1.5390895594100267e-09; }
    consteval double expressionFrequency(Isotope<170>) noexcept { return 1.026059992426501e-09; }
    consteval double expressionFrequency(Isotope<171>) noexcept { return 6.840399735422548e-10; }
    consteval double expressionFrequency(Isotope<172>) noexcept { return 1.026059929155679e-09; }
    consteval double expressionFrequency(Isotope<173>) noexcept { return 8.550500955701031e-10; }
    consteval double expressionFrequency(Isotope<174>) noexcept { return 8.550499374263083e-10; }
    consteval double expressionFrequency(Isotope<175>) noexcept { return 8.550497916649192e-10; }
    consteval double expressionFrequency(Isotope<176>) noexcept { return 8.550500955701031e-10; }
    consteval double expressionFrequency(Isotope<177>) noexcept { return 2.206027561103357e-08; }
    consteval double expressionFrequency(Isotope<178>) noexcept { return 1.026059866840272e-09; }
    consteval double expressionFrequency(Isotope<179>) noexcept { return 6.84039911226848e-10; }
    consteval double expressionFrequency(Isotope<180>) noexcept { return 8.550497986703092e-10; }
    consteval double expressionFrequency(Isotope<181>) noexcept { return 1.539090077092898e-09; }
    consteval double expressionFrequency(Isotope<182>) noexcept { return 1.026059929155679e-09; }
    consteval double expressionFrequency(Isotope<183>) noexcept { return 5.130300580509495e-10; }
    consteval double expressionFrequency(Isotope<184>) noexcept { return 1.0260599998381425e-09; }
    consteval double expressionFrequency(Isotope<185>) noexcept { return 8.550499513489667e-10; }
    consteval double expressionFrequency(Isotope<186>) noexcept { return 5.643328885711218e-09; }
    consteval double expressionFrequency(Isotope<187>) noexcept { return 1.0260601265382283e-09; }
    consteval double expressionFrequency(Isotope<188>) noexcept { return 8.037475363809575e-09; }
    consteval double expressionFrequency(Isotope<189>) noexcept { return 6.84039911226848e-10; }
    consteval double expressionFrequency(Isotope<190>) noexcept { return 5.130301370019685e-10; }
    consteval double expressionFrequency(Isotope<191>) noexcept { return 1.5390890189593477e-09; }
    consteval double expressionFrequency(Isotope<192>) noexcept { return 6.840399435322592e-10; }
    consteval double expressionFrequency(Isotope<193>) noexcept { return 1.7101019039135672e-09; }
    consteval double expressionFrequency(Isotope<194>) noexcept { return 1.5390898825821535e-09; }
    consteval double expressionFrequency(Isotope<195>) noexcept { return 3.078181742949178e-09; }
    consteval double expressionFrequency(Isotope<196>) noexcept { return 8.550499204527973e-10; }
    consteval double expressionFrequency(Isotope<197>) noexcept { return 4.104238373321323e-09; }
    consteval double expressionFrequency(Isotope<198>) noexcept { return 2.0521199055274886e-09; }
    consteval double expressionFrequency(Isotope<199>) noexcept { return 1.3680797203233845e-09; }
    consteval double expressionFrequency(Isotope<200>) noexcept { return 1.3680806136661872e-09; }
    consteval double expressionFrequency(Isotope<201>) noexcept { return 3.5912089792006946e-09; }
    consteval double expressionFrequency(Isotope<202>) noexcept { return 1.0260608320096806e-09; }
    consteval double expressionFrequency(Isotope<203>) noexcept { return 4.9592901637446496e-09; }
    consteval double expressionFrequency(Isotope<204>) noexcept { return 1.1970697811076243e-09; }
    consteval double expressionFrequency(Isotope<205>) noexcept { return 1.0260608180311017e-09; }
    consteval double expressionFrequency(Isotope<206>) noexcept { return 6.840407424133684e-10; }
    consteval double expressionFrequency(Isotope<207>) noexcept { return 8.550494152456521e-10; }
    consteval double expressionFrequency(Isotope<208>) noexcept { return 1.0260592560600596e-09; }
    consteval double expressionFrequency(Isotope<209>) noexcept { return 1.3680796509172029e-09; }
    consteval double expressionFrequency(Isotope<210>) noexcept { return 1.7101000789968155e-09; }
    consteval double expressionFrequency(Isotope<211>) noexcept { return 1.983716034711261e-08; }
    consteval double expressionFrequency(Isotope<212>) noexcept { return 1.5390896023473345e-09; }
    consteval double expressionFrequency(Isotope<213>) noexcept { return 1.7101007451501643e-09; }
    consteval double expressionFrequency(Isotope<214>) noexcept { return 5.1302973703596e-10; }
    consteval double expressionFrequency(Isotope<215>) noexcept { return 1.7101009102462676e-10; }
    consteval double expressionFrequency(Isotope<216>) noexcept { return 3.420202062841391e-10; }
    consteval double expressionFrequency(Isotope<217>) noexcept { return 5.13030203518817e-10; }
    consteval double expressionFrequency(Isotope<218>) noexcept { return 8.550499574537838e-10; }
    consteval double expressionFrequency(Isotope<219>) noexcept { return 1.71009977806712e-10; }
    consteval double expressionFrequency(Isotope<220>) noexcept { return 8.550499306169814e-10; }
    consteval double expressionFrequency(Isotope<221>) noexcept { return 1.7101004848917573e-10; }
    consteval double expressionFrequency(Isotope<222>) noexcept { return 5.130301089534575e-10; }
    consteval double expressionFrequency(Isotope<223>) noexcept { return 5.130295126898212e-10; }
    consteval double expressionFrequency(Isotope<224>) noexcept { return 5.130300466380508e-10; }
    consteval double expressionFrequency(Isotope<225>) noexcept { return 5.130305035149808e-09; }
    consteval double expressionFrequency(Isotope<226>) noexcept { return 6.840400843953236e-10; }
    consteval double expressionFrequency(Isotope<227>) noexcept { return 1.5390903646535457e-09; }
    consteval double expressionFrequency(Isotope<228>) noexcept { return 1.0260595902152e-09; }
    consteval double expressionFrequency(Isotope<229>) noexcept { return 3.420200337924868e-10; }
    consteval double expressionFrequency(Isotope<230>) noexcept { return 8.550500022514747e-10; }
    consteval double expressionFrequency(Isotope<231>) noexcept { return 6.84039911226848e-10; }
    consteval double expressionFrequency(Isotope<232>) noexcept { return 6.840400486796483e-10; }
    consteval double expressionFrequency(Isotope<233>) noexcept { return 8.550496272704871e-10; }
    consteval double expressionFrequency(Isotope<234>) noexcept { return 1.0260599385077532e-09; }
    consteval double expressionFrequency(Isotope<235>) noexcept { return 1.197069477205724e-09; }
    consteval double expressionFrequency(Isotope<236>) noexcept { return 6.840397401697642e-10; }
    consteval double expressionFrequency(Isotope<237>) noexcept { return 1.368079764614527e-09; }
    consteval double expressionFrequency(Isotope<238>) noexcept { return 1.7101000975311516e-09; }
    consteval double expressionFrequency(Isotope<239>) noexcept { return 1.8811093074129653e-09; }
    consteval double expressionFrequency(Isotope<240>) noexcept { return 8.5504988903356e-10; }
    consteval double expressionFrequency(Isotope<241>) noexcept { return 1.5390894079875917e-09; }
    consteval double expressionFrequency(Isotope<242>) noexcept { return 5.13029933420136e-10; }
    consteval double expressionFrequency(Isotope<243>) noexcept { return 5.13029933420136e-10; }
    consteval double expressionFrequency(Isotope<244>) noexcept { return 8.550500022514747e-10; }
    consteval double expressionFrequency(Isotope<245>) noexcept { return 6.84039911226848e-10; }
    consteval double expressionFrequency(Isotope<246>) noexcept { return 8.5504988903356e-10; }
    consteval double expressionFrequency(Isotope<247>) noexcept { return 3.42019955613424e-10; }
    consteval double expressionFrequency(Isotope<248>) noexcept { return 5.13029933420136e-10; }
    consteval double expressionFrequency(Isotope<249>) noexcept { return 1.8811105051128452e-09; }
    consteval double expressionFrequency(Isotope<250>) noexcept { return 6.84039911226848e-10; }
    consteval double expressionFrequency(Isotope<251>) noexcept { return 3.42019955613424e-10; }
    consteval double expressionFrequency(Isotope<252>) noexcept { return 8.55050337543917e-10; }
    consteval double expressionFrequency(Isotope<253>) noexcept { return 1.0260599640488474e-09; }
    consteval double expressionFrequency(Isotope<254>) noexcept { return 5.13029933420136e-10; }
    consteval double expressionFrequency(Isotope<255>) noexcept { return 1.71009977806712e-10; }
    consteval double expressionFrequency(Isotope<256>) noexcept { return 0.0; }
    consteval double expressionFrequency(Isotope<257>) noexcept { return 3.42019955613424e-10; }
    consteval double expressionFrequency(Isotope<258>) noexcept { return 1.71009977806712e-10; }
    consteval double expressionFrequency(Isotope<259>) noexcept { return 0.0; }
    consteval double expressionFrequency(Degree<0>) noexcept { return 0.00039172405592891144; }
    consteval double expressionFrequency(Degree<1>) noexcept { return 0.5201466586307003; }
    consteval double expressionFrequency(Degree<2>) noexcept { return 0.04864862411087598; }
    consteval double expressionFrequency(Degree<3>) noexcept { return 0.2676745431637454; }
    consteval double expressionFrequency(Degree<4>) noexcept { return 0.1631255747936151; }
    consteval double expressionFrequency(Degree<5>) noexcept { return 6.696240280797756e-06; }
    consteval double expressionFrequency(Degree<6>) noexcept { return 6.161830782656515e-06; }
    consteval double expressionFrequency(Degree<7>) noexcept { return 2.6335533544154077e-08; }
    consteval double expressionFrequency(Degree<8>) noexcept { return 8.550504835318602e-10; }
    consteval double expressionFrequency(Degree<9>) noexcept { return 0.0; }
    consteval double expressionFrequency(Valence<0>) noexcept { return 0.00039172405592891144; }
    consteval double expressionFrequency(Valence<1>) noexcept { return 0.49070386818883854; }
    consteval double expressionFrequency(Valence<2>) noexcept { return 0.0594899825012628; }
    consteval double expressionFrequency(Valence<3>) noexcept { return 0.05172710293294454; }
    consteval double expressionFrequency(Valence<4>) noexcept { return 0.39510366353926923; }
    consteval double expressionFrequency(Valence<5>) noexcept { return 0.0003466442016566324; }
    consteval double expressionFrequency(Valence<6>) noexcept { return 0.002231680314681798; }
    consteval double expressionFrequency(Valence<7>) noexcept { return 5.3187523750552235e-06; }
    consteval double expressionFrequency(Valence<8>) noexcept { return 2.3770376050354293e-08; }
    consteval double expressionFrequency(Valence<9>) noexcept { return 3.4201994106459537e-10; }
    consteval double expressionFrequency(Connectivity<0>) noexcept { return 0.00039172405592891144; }
    consteval double expressionFrequency(Connectivity<1>) noexcept { return 0.5201466586307003; }
    consteval double expressionFrequency(Connectivity<2>) noexcept { return 0.04864862411087598; }
    consteval double expressionFrequency(Connectivity<3>) noexcept { return 0.2676745431637454; }
    consteval double expressionFrequency(Connectivity<4>) noexcept { return 0.1631255747936151; }
    consteval double expressionFrequency(Connectivity<5>) noexcept { return 6.696240280797756e-06; }
    consteval double expressionFrequency(Connectivity<6>) noexcept { return 6.161830782656515e-06; }
    consteval double expressionFrequency(Connectivity<7>) noexcept { return 2.6335533544154077e-08; }
    consteval double expressionFrequency(Connectivity<8>) noexcept { return 8.550504835318602e-10; }
    consteval double expressionFrequency(Connectivity<9>) noexcept { return 0.0; }
    consteval double expressionFrequency(TotalH<0>) noexcept { return 0.697546654786409; }
    consteval double expressionFrequency(TotalH<1>) noexcept { return 0.1704494158883208; }
    consteval double expressionFrequency(TotalH<2>) noexcept { return 0.09158444925107122; }
    consteval double expressionFrequency(TotalH<3>) noexcept { return 0.040355017069505485; }
    consteval double expressionFrequency(TotalH<4>) noexcept { return 6.442329120397767e-05; }
    consteval double expressionFrequency(TotalH<5>) noexcept { return 3.7622201359778625e-09; }
    consteval double expressionFrequency(TotalH<6>) noexcept { return 3.249188358317651e-09; }
    consteval double expressionFrequency(TotalH<7>) noexcept { return 1.881107913193674e-09; }
    consteval double expressionFrequency(TotalH<8>) noexcept { return 0.0; }
    consteval double expressionFrequency(TotalH<9>) noexcept { return 0.0; }
    consteval double expressionFrequency(ImplicitH<0>) noexcept { return 0.697546654786409; }
    consteval double expressionFrequency(ImplicitH<1>) noexcept { return 0.1704494158883208; }
    consteval double expressionFrequency(ImplicitH<2>) noexcept { return 0.09158444925107122; }
    consteval double expressionFrequency(ImplicitH<3>) noexcept { return 0.040355017069505485; }
    consteval double expressionFrequency(ImplicitH<4>) noexcept { return 6.442329120397767e-05; }
    consteval double expressionFrequency(ImplicitH<5>) noexcept { return 3.7622201359778625e-09; }
    consteval double expressionFrequency(ImplicitH<6>) noexcept { return 3.249188358317651e-09; }
    consteval double expressionFrequency(ImplicitH<7>) noexcept { return 1.881107913193674e-09; }
    consteval double expressionFrequency(ImplicitH<8>) noexcept { return 0.0; }
    consteval double expressionFrequency(ImplicitH<9>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingCount<0>) noexcept { return 0.6947270642281514; }
    consteval double expressionFrequency(RingCount<1>) noexcept { return 0.27944265779913896; }
    consteval double expressionFrequency(RingCount<2>) noexcept { return 0.02509977548166723; }
    consteval double expressionFrequency(RingCount<3>) noexcept { return 0.0006879928036211851; }
    consteval double expressionFrequency(RingCount<4>) noexcept { return 4.022735976134487e-05; }
    consteval double expressionFrequency(RingCount<5>) noexcept { return 1.9760207836528006e-06; }
    consteval double expressionFrequency(RingCount<6>) noexcept { return 2.175246366365939e-07; }
    consteval double expressionFrequency(RingCount<7>) noexcept { return 5.472319984109799e-08; }
    consteval double expressionFrequency(RingCount<8>) noexcept { return 1.983715446525459e-08; }
    consteval double expressionFrequency(RingCount<9>) noexcept { return 6.327371224155653e-09; }
    consteval double expressionFrequency(RingSize<0>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingSize<1>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingSize<2>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingSize<3>) noexcept { return 0.003456892015976178; }
    consteval double expressionFrequency(RingSize<4>) noexcept { return 0.002574508620671572; }
    consteval double expressionFrequency(RingSize<5>) noexcept { return 0.07137185992930796; }
    consteval double expressionFrequency(RingSize<6>) noexcept { return 0.23395455007659788; }
    consteval double expressionFrequency(RingSize<7>) noexcept { return 0.004644510022728279; }
    consteval double expressionFrequency(RingSize<8>) noexcept { return 0.00069611069185377; }
    consteval double expressionFrequency(RingSize<9>) noexcept { return 0.0002172745150120276; }
    consteval double expressionFrequency(RingConnectivity<0>) noexcept { return 0.6946352665958881; }
    consteval double expressionFrequency(RingConnectivity<1>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingConnectivity<2>) noexcept { return 0.2805731166955411; }
    consteval double expressionFrequency(RingConnectivity<3>) noexcept { return 0.0241919006083607; }
    consteval double expressionFrequency(RingConnectivity<4>) noexcept { return 0.0005995780697728736; }
    consteval double expressionFrequency(RingConnectivity<5>) noexcept { return 8.34528467584299e-08; }
    consteval double expressionFrequency(RingConnectivity<6>) noexcept { return 6.412877267236231e-08; }
    consteval double expressionFrequency(RingConnectivity<7>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingConnectivity<8>) noexcept { return 0.0; }
    consteval double expressionFrequency(RingConnectivity<9>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-15>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-14>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-13>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-12>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-11>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-10>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-9>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-8>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-7>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-6>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-5>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<-4>) noexcept { return 1.7100994873650872e-09; }
    consteval double expressionFrequency(Charge<-3>) noexcept { return 6.994307315866859e-08; }
    consteval double expressionFrequency(Charge<-2>) noexcept { return 7.426452439252166e-06; }
    consteval double expressionFrequency(Charge<-1>) noexcept { return 0.001584548612321732; }
    consteval double expressionFrequency(Charge<0>) noexcept { return 0.9968068696774753; }
    consteval double expressionFrequency(Charge<1>) noexcept { return 0.0015222022197550924; }
    consteval double expressionFrequency(Charge<2>) noexcept { return 5.628382032097709e-05; }
    consteval double expressionFrequency(Charge<3>) noexcept { return 1.4406563252132019e-05; }
    consteval double expressionFrequency(Charge<4>) noexcept { return 7.357191452468266e-06; }
    consteval double expressionFrequency(Charge<5>) noexcept { return 4.605300823946001e-07; }
    consteval double expressionFrequency(Charge<6>) noexcept { return 3.300490922961394e-07; }
    consteval double expressionFrequency(Charge<7>) noexcept { return 2.4796453233565955e-08; }
    consteval double expressionFrequency(Charge<8>) noexcept { return 2.5309464280626407e-08; }
    consteval double expressionFrequency(Charge<9>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<10>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<11>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<12>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<13>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<14>) noexcept { return 0.0; }
    consteval double expressionFrequency(Charge<15>) noexcept { return 0.0; }
    consteval double expressionFrequency(AromaticAtom<5>) noexcept { return 2.470376461502374e-05; }
    consteval double expressionFrequency(AromaticAtom<6>) noexcept { return 0.1974038913919688; }
    consteval double expressionFrequency(AromaticAtom<7>) noexcept { return 0.0193999658042352; }
    consteval double expressionFrequency(AromaticAtom<8>) noexcept { return 0.0017222683454965564; }
    consteval double expressionFrequency(AromaticAtom<15>) noexcept { return 4.275933532580179e-06; }
    consteval double expressionFrequency(AromaticAtom<16>) noexcept { return 0.002299238152721317; }
    consteval double expressionFrequency(AromaticAtom<33>) noexcept { return 1.0192196392820635e-07; }
    consteval double expressionFrequency(AromaticAtom<34>) noexcept { return 4.755787892217665e-06; }
    consteval double expressionFrequency(AnyAtom) noexcept { return 1.0; }
    consteval double expressionFrequency(AnyAliphatic) noexcept { return 0.7791408016780832; }
    consteval double expressionFrequency(AnyAromatic) noexcept { return 0.22085919979187205; }
    consteval double expressionFrequency(Cyclic) noexcept { return 0.3053647342517278; }
    consteval double expressionFrequency(Acyclic) noexcept { return 0.6946352665958881; }

    // Bond Primitives

    consteval double expressionFrequency(ImplicitBond) noexcept { return 0.9635761238622592; }
    consteval double expressionFrequency(AnyBond) noexcept { return 1.0; }
    consteval double expressionFrequency(AromaticBond) noexcept { return 0.21854685393206996; }
    consteval double expressionFrequency(RingBond) noexcept { return 0.30685414335650296; }
    consteval double expressionFrequency(BondOrder<1>) noexcept { return 0.7450293257198629; }
    consteval double expressionFrequency(BondOrder<2>) noexcept { return 0.03499309764249357; }
    consteval double expressionFrequency(BondOrder<3>) noexcept { return 0.001430744492880069; }
    consteval double expressionFrequency(BondOrder<4>) noexcept { return 0.0; }

} // namespace Kitimar::CTSmarts

#include <ctll/list.hpp>

#include <algorithm>
#include <functional>

namespace Kitimar::CTSmarts {

    //
    // Sort ctll::list directly
    //

    template<typename Project = std::identity, typename Compare = std::less<void>>
    consteval auto selectLast(auto input)
    {
        if constexpr (ctll::size(input) == 1)
            return ctll::front(input);
        else {
            auto [a, tail] = ctll::pop_and_get_front(input);
            auto b = selectLast<Project, Compare>(tail);
            if constexpr (Compare{}(Project{}(a), Project{}(b)))
                return b;
            else
                return a;
        }
    }

    template<typename Project = std::identity, typename Compare = std::less<>, typename Output = ctll::empty_list>
    consteval auto selectionSort(auto input, Output output = {})
    {
        if constexpr (!ctll::empty(input)) {
            auto last = selectLast<Project, Compare>(input);
            return selectionSort<Project, Compare>(ctll::remove_item(last, input), ctll::push_front(last, output));
        } else
            return output;
    }

    //
    // Sort ctll::list indirectly
    //
    // ctll::list<Expr...> -> transform -> std::array< [ Index, Project(Expr) ] > -> std::sort -> transform -> ctll::list<Expr...>

    /*
    template<typename Project, typename T, int N>
    constexpr auto makeSortableHelper(auto input)
    {
        if constexpr (!ctll::empty(input)) {
            auto [head, tail] = ctll::pop_and_get_front(input);
            auto sortable = makeSortableHelper<Project, T, N>(tail);
            auto index = N - ctll::size(input);
            sortable[index] = std::make_pair(index, Project{}(head));
            return sortable;
        } else
            return std::array<std::pair<int, T>, N>{};
    }


    template<typename Project>
    constexpr auto makeSortable(auto input)
    {
        using T = decltype(Project{}(ctll::front(input)));
        constexpr auto N = ctll::size(input);
        return makeSortableHelper<Project, T, N>(input);
    }


    template<typename SortableT, typename Compare>
    constexpr auto sortSortable()
    {
        //return SortableT::data;

        auto copy = SortableT::data;
        std::ranges::sort(copy, [] (const auto &a, const auto &b) {
            return !Compare{}(a.second, b.second);
        });
        return copy;

    }

    template<typename Input, typename Project>
    struct Sortable
    {
        static constexpr inline auto data = makeSortable<Project>(Input{});

        consteval Sortable() noexcept {}
        consteval Sortable(Input, Project) noexcept {}
    };

    template<typename SortableT, typename Compare>
    struct SortedSortable
    {
        static constexpr inline auto data = sortSortable<SortableT, Compare>();

        consteval SortedSortable() noexcept {}
        consteval SortedSortable(SortableT, Compare) noexcept {}
    };

    template<typename SortedT, int I = 0>
    constexpr auto makeSorted(auto input)
    {
        if constexpr (I == ctll::size(input))
            return ctll::empty_list{};
        else {
            constexpr auto sortableIndex = ctll::size(input) - I - 1;
            constexpr auto inputIndex = SortedT::data[sortableIndex].first;
            return ctll::push_front(get<inputIndex>(input), makeSorted<SortedT, I + 1>(input));
        }
    }

    template<typename Input, typename SortableT>
    struct Sorted
    {
        static constexpr inline auto data = makeSorted<SortableT>(Input{});

        consteval Sorted() noexcept {}
        consteval Sorted(Input, SortableT) noexcept {}
    };

    template<typename Project = std::identity, typename Compare = std::less<>>
    constexpr auto stdSort(auto input)
    {
        //auto sortable = makeSortable<Project>(input);
        //sortSortable<Compare>(sortable);
        auto sortable = Sortable{input, Project{}};
        auto sortedSortable = SortedSortable{sortable, Compare{}};
        auto sorted = Sorted{input, sortedSortable};
        return sorted.data;
    }
    */

    template<typename Project = std::identity, typename Compare = std::less<>>
    consteval auto ctllSort(auto input)
    {
        return selectionSort<Project, Compare>(input);
        //return stdSort<Project, Compare>(input);
    }

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    struct ProjExprFrequency;

    template<typename Expr>
    consteval double expressionFrequency(Not<Expr> op) noexcept
    {
        return 1 - expressionFrequency(op.expr);
    }

    template<typename ...Expr>
    consteval double expressionFrequency(And<Expr...> op) noexcept
    {        
        return expressionFrequency(selectLast<ProjExprFrequency, std::greater<>>(op.expr));
    }

    template<typename ...Expr>
    consteval double expressionFrequency(Or<Expr...> op) noexcept
    {        
        return expressionFrequency(selectLast<ProjExprFrequency>(op.expr));
    }

    consteval double expressionFrequency(...) noexcept { return 0; }

    struct ProjExprFrequency
    {
        consteval auto operator()(auto expr)
        {
            return expressionFrequency(expr);
        }
    };

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    struct BestCase {};
    struct WorstCase {};

    // Not

    template<typename Goal = BestCase, typename Expr>
    consteval auto optimizeExpression(Not<Expr> op)
    {
        using NotGoal = std::conditional_t<std::is_same_v<Goal, BestCase>, WorstCase, BestCase>;
        return Not(optimizeExpression<NotGoal>(op.expr));
    }

    // And

    template<typename Goal = BestCase, typename ...Expr>
    consteval auto optimizeExpression(And<Expr...> op)
    {        
        using Compare = std::conditional_t<std::is_same_v<Goal, BestCase>, std::less<>, std::greater<>>;
        return And(ctllSort<ProjExprFrequency, Compare>(optimizeExpressions<Goal>(op.expr)));
    }

    // Or

    template<typename Goal = BestCase, typename ...Expr>
    consteval auto optimizeExpression(Or<Expr...> op)
    {
        using Compare = std::conditional_t<std::is_same_v<Goal, BestCase>, std::greater<>, std::less<>>;
        return Or(ctllSort<ProjExprFrequency, Compare>(optimizeExpressions<Goal>(op.expr)));
    }

    // Primitive

    template<typename Goal = BestCase>
    consteval auto optimizeExpression(auto expr)
    {
        return expr;
    }

    // ctll::list<Expr>

    template<typename Goal = BestCase, typename ...Expr>
    consteval auto optimizeExpressions(ctll::list<Expr...> expressions) noexcept
    {
        if constexpr (!ctll::empty(expressions)) {
            auto [expr, tail] = ctll::pop_and_get_front(expressions);
            return ctll::push_front(optimizeExpression<Goal>(expr), optimizeExpressions<Goal>(tail));
        } else
            return ctll::empty_list{};
    }

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    struct ProjAtomFrequency
    {        
        consteval auto operator()(auto atom)
        {
            return expressionFrequency(atom.expr);
        }
    };

    template<typename SmartsT>
    consteval auto makeAtomFrequency(auto atoms) noexcept
    {
        if constexpr (ctll::empty(atoms))
            return std::array<double, SmartsT::numAtoms>{};
        else {
            auto [atom, tail] = ctll::pop_and_get_front(atoms);
            auto freq = makeAtomFrequency<SmartsT>(tail);
            freq[atom.index] = expressionFrequency(atom.expr);
            return freq;
        }
    }

    template<typename SmartsT>
    struct AtomFrequency
    {
        static constexpr inline auto data = makeAtomFrequency<SmartsT>(SmartsT::atoms);

        consteval AtomFrequency() noexcept {}
        consteval AtomFrequency(SmartsT) noexcept {}
    };

} // namespace Kitimar::CTSmarts

#include <algorithm>
#include <type_traits>

namespace Kitimar::CTSmarts {

    template<typename SmartsT, typename IncidentListT, typename VertexFrequencyT, bool SeedBond>
    consteval auto makeOptimizeIncidentList()
    {
        auto incident = IncidentListT::data;

        for (auto i = 0; i < static_cast<int>(SmartsT::numAtoms); ++i) {
            auto offset = i * IncidentListT::stride;
            auto end = offset + IncidentListT::degrees.data[i];
            if (SeedBond)
                ++offset;
            std::ranges::sort(std::begin(incident) + offset, std::begin(incident) + end, [i] (auto index1, auto index2) {
                auto edge1 = IncidentListT::edges.data[index1];
                auto edge2 = IncidentListT::edges.data[index2];
                auto target1 = i == edge1.source ? edge1.target : edge1.source;
                auto target2 = i == edge2.source ? edge2.target : edge2.source;
                return VertexFrequencyT::data[target1] < VertexFrequencyT::data[target2];
            });
        }

        return incident;
    }

    template<typename SmartsT, typename IncidentListT, typename VertexFrequencyT, bool SeedBond>
    struct OptimizeIncidentList
    {
        // store adjacent (or incident) bond indices for each vertex
        static constexpr inline auto data = makeOptimizeIncidentList<SmartsT, IncidentListT, VertexFrequencyT, SeedBond>();
        static constexpr inline auto edges = IncidentListT::edges;
        static constexpr inline auto degrees = IncidentListT::degrees;
        static constexpr inline auto stride = IncidentListT::stride;

        static consteval auto get(int AtomIndex, int AdjIndex)
        {
            return data[stride * AtomIndex + AdjIndex];
        }

        consteval OptimizeIncidentList() noexcept {}
        consteval OptimizeIncidentList(SmartsT, IncidentListT, VertexFrequencyT, std::integral_constant<bool, SeedBond>) noexcept {}
    };

}

#include <ctll/list.hpp>
#include <ctll/grammars.hpp>

namespace Kitimar::CTSmarts {

    struct SmartsGrammar
    {
        template<bool Component = false, typename Recursive = ctll::empty_list, int Branch = 0>
        struct Parenthesis
        {
            static constexpr inline auto component = Component;
            static constexpr inline auto recursive = Recursive{};
            static constexpr inline auto branch = Branch;
        };

        template<typename P>
        using ParenthesisPushComponent = Parenthesis<true, decltype(P::recursive), P::branch>;

        template<typename P>
        using ParenthesisPopComponent = Parenthesis<false, decltype(P::recursive), P::branch>;

        template<typename P>
        using ParenthesisPushRecursive = Parenthesis<P::component, decltype(ctll::push_front(Number<P::branch>{}, P::recursive)), 0>;

        template<typename P>
        using ParenthesisPopRecursive = Parenthesis<P::component, decltype(ctll::pop_front(P::recursive)), ctll::front(P::recursive).value>;

        template<typename P>
        using ParenthesisPushBranch = Parenthesis<P::component, decltype(P::recursive), P::branch + 1>;

        template<typename P>
        using ParenthesisPopBranch = Parenthesis<P::component, decltype(P::recursive), P::branch - 1>;

        //
        // Symbols
        //

        template<typename P> struct atom {};
        template<typename P> struct atom_B {};
        template<typename P> struct atom_C {};
        template<typename P> struct atom_expr {}; // FIXME: rename to bracket_atom
        template<typename P> struct atom_expr2 {}; // FIXME: rename to atom_expr
        template<typename P> struct atom_exprA {}; // atom_expr_A
        template<typename P> struct atom_exprB {};
        template<typename P> struct atom_exprC {};
        template<typename P> struct atom_exprD {};
        template<typename P> struct atom_exprE {};
        template<typename P> struct atom_exprF {};
        template<typename P> struct atom_exprG {};
        template<typename P> struct atom_exprH {};
        template<typename P> struct atom_exprI {};
        template<typename P> struct atom_exprK {};
        template<typename P> struct atom_exprL {};
        template<typename P> struct atom_exprM {};
        template<typename P> struct atom_exprN {};
        template<typename P> struct atom_exprO {};
        template<typename P> struct atom_exprP {};
        template<typename P> struct atom_exprR {};
        template<typename P> struct atom_exprS {};
        template<typename P> struct atom_exprT {};
        template<typename P> struct atom_exprX {};
        template<typename P> struct atom_exprY {};
        template<typename P> struct atom_exprZ {};
        template<typename P> struct atom_expr_a {};
        template<typename P> struct atom_expr_s {};
        template<typename P> struct atom_expr_isotope {};
        template<typename P> struct atom_expr_isotope2 {};
        template<typename P> struct atom_expr_element {};
        template<typename P> struct atom_expr_element2 {};
        template<typename P> struct atom_expr_degree {};
        template<typename P> struct atom_expr_valence {};
        template<typename P> struct atom_expr_valence2 {};
        template<typename P> struct atom_expr_connectivity {};
        template<typename P> struct atom_expr_total_h {};
        template<typename P> struct atom_expr_impl_h {};
        template<typename P> struct atom_expr_ring_count {};
        template<typename P> struct atom_expr_ring_size {};
        template<typename P> struct atom_expr_ring_size2 {};
        template<typename P> struct atom_expr_ring_connectivity {};
        template<typename P> struct atom_expr_ring_connectivity2 {};
        template<typename P> struct atom_expr_neg_charge {};
        template<typename P> struct atom_expr_neg_charge2 {};
        template<typename P> struct atom_expr_pos_charge {};
        template<typename P> struct atom_expr_pos_charge2 {};
        template<typename P> struct atom_expr_chiral {};
        template<typename P> struct atom_expr_class {};
        template<typename P> struct atom_expr_class2 {};
        template<typename P> struct atom_expr_recursive {};
        template<typename P> struct bond_expr {};
        template<typename P> struct bond_expr2 {};
        template<typename P> struct chain {};
        template<typename P> struct chain_up_down {};
        template<typename P> struct ring_bond {};
        template<typename P> struct ring_bond2 {};

        //
        // Actions
        //

        // Create atom AST elements
        struct make_any_atom : ctll::action {}; // '*'
        struct make_any_aliphatic : ctll::action {}; // 'A'
        struct make_any_aromatic : ctll::action {}; // 'a'
        struct make_aliphatic : ctll::action {}; // single letter symbols
        struct make_aromatic : ctll::action {}; // single letter symbols

        struct make_isotope : ctll::action {};
        struct make_element : ctll::action {};
        struct make_degree : ctll::action {};
        struct make_valence : ctll::action {};
        struct make_connectivity : ctll::action {};
        struct make_total_h : ctll::action {};
        struct make_has_impl_h : ctll::action {};
        struct make_impl_h : ctll::action {};
        struct make_cyclic : ctll::action {};
        struct make_acyclic : ctll::action {};
        struct make_ring_count : ctll::action {};
        struct make_ring_size : ctll::action {};
        struct make_ring_connectivity : ctll::action {};
        struct start_charge : ctll::action {};
        struct increment_charge : ctll::action {};
        struct decrement_charge : ctll::action {};
        struct make_charge : ctll::action {};
        struct make_chiral : ctll::action {};
        struct make_class : ctll::action {};

        struct make_bond_primitive : ctll::action {};
        //struct make_up_or_down_bond : ctll::action {};

        // Create bond AST elements
        struct next_atom : ctll::action {};

        // Components
        struct push_component : ctll::action {};
        struct pop_component : ctll::action {};

        // Branches
        struct push_prev : ctll::action {};
        struct pop_prev : ctll::action {};
        struct reset_prev : ctll::action {};
        struct set_bond_type : ctll::action {};
        struct handle_ring_bond : ctll::action {};

        // Recursive
        struct push_recursive : ctll::action {};
        struct pop_recursive : ctll::action {};

        // Operators
        /*
        struct toggle_not : ctll::action {};
        struct reset_not : ctll::action {};
        struct set_and_high : ctll::action {};
        struct set_or : ctll::action {};
        struct set_and_low : ctll::action {};
        struct set_no_op : ctll::action {};
        */

        struct make_atom_not : ctll::action {};
        struct make_atom_and_high : ctll::action {};
        struct make_atom_or : ctll::action {};
        struct make_atom_and_low : ctll::action {};

        struct make_bond_not : ctll::action {};
        struct make_bond_and_high : ctll::action {};
        struct make_bond_or : ctll::action {};
        struct make_bond_and_low : ctll::action {};

        // Helpers
        struct push_char : ctll::action {};
        struct pop_char : ctll::action {};

        struct start_number : ctll::action {};
        struct push_number : ctll::action {};

        struct error_empty_bracket : ctll::action {}; // '[]'

        using _start = atom<Parenthesis<>>;

        //
        // Organic atoms
        //

        // aliphatic atom: 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'N', 'O', 'P', 'S', 'F', 'I'>) -> ctll::push<ctll::anything, make_aliphatic, next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'B'>) -> ctll::push<ctll::anything, push_char, atom_B<P>>;
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'C'>) -> ctll::push<ctll::anything, push_char, atom_C<P>>;
        template<typename P> static constexpr auto rule(atom_B<P>, ctll::term<'r'>) -> ctll::push<ctll::anything, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_B<P>, ctll::neg_set<'r'>) -> ctll::push<pop_char, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_B<P>, ctll::epsilon) -> ctll::push<pop_char, make_aliphatic, next_atom>;
        template<typename P> static constexpr auto rule(atom_C<P>, ctll::term<'l'>) -> ctll::push<ctll::anything, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_C<P>, ctll::neg_set<'l'>) -> ctll::push<pop_char, make_aliphatic,  next_atom, chain<P>>;
        template<typename P> static constexpr auto rule(atom_C<P>, ctll::epsilon) -> ctll::push<pop_char, make_aliphatic, next_atom>;
        // aromatic atom: 'b' | 'c' | 'n' | 'o' | 's' | 'p'
        template<typename P> static constexpr auto rule(atom<P>, ctll::set<'b', 'c','n','o','p','s'>) -> ctll::push<ctll::anything, make_aromatic, next_atom, chain<P>>;
        // any atom: '*'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'*'>) -> ctll::push<ctll::anything, make_any_atom, next_atom, chain<P>>;
        // any aromatic: 'a'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'a'>) -> ctll::push<ctll::anything, make_any_aromatic, next_atom, chain<P>>;
        // any aliphatic: 'A'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'A'>) -> ctll::push<ctll::anything, make_any_aliphatic, next_atom, chain<P>>;
        // bracket atom: '[' atom_expression+ ']'
        template<typename P> static constexpr auto rule(atom<P>, ctll::term<'['>) -> ctll::push<ctll::anything, atom_expr<P>>;

        template<typename P> static constexpr auto rule(atom<P>, ctll::neg_set<'B','C','N','O','P','S','F','I','b','c','n','o','p','s','*','a','A','['>) -> ctll::reject;

        //
        // Atom expressions (i.e. [ ... ] )
        //

        using digit_chars = ctll::range<'0', '9'>;
        using not_digit_chars = ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>;

        // check for [], [!] or [...!]
        template<typename P> static constexpr auto rule(atom_expr<P>, ctll::set<']'>) -> ctll::push<error_empty_bracket, ctll::reject>;
        template<typename P> static constexpr auto rule(atom_expr<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_atom_not, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr<P>, ctll::neg_set<']', '!'>) -> ctll::push<atom_expr2<P>>;

        // operations + end ]
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_atom_not, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'&'>) -> ctll::push<ctll::anything, make_atom_and_high, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<','>) -> ctll::push<ctll::anything, make_atom_or, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<';'>) -> ctll::push<ctll::anything, make_atom_and_low, atom_expr<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<']'>) -> ctll::push<ctll::anything, next_atom, chain<P>>;

        // isotope: NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_isotope<P>>;
        template<typename P> static constexpr auto rule(atom_expr_isotope<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_isotope<P>>;
        template<typename P> static constexpr auto rule(atom_expr_isotope<P>, not_digit_chars) -> ctll::push<make_isotope, atom_expr2<P>>;

        // element: '#' | '#' NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'#'>) -> ctll::push<ctll::anything, atom_expr_element<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_element2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element<P>, not_digit_chars) -> ctll::push<make_element, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_element2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_element2<P>, not_digit_chars) -> ctll::push<make_element, atom_expr2<P>>;

        // valence: 'v' | 'v' NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'v'>) -> ctll::push<ctll::anything, atom_expr_valence<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_valence2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence<P>, not_digit_chars) -> ctll::push<make_valence, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_valence2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_valence2<P>, not_digit_chars) -> ctll::push<make_valence, atom_expr2<P>>;

        // implicit hydrogens: 'h' | 'h' DIGIT
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'h'>) -> ctll::push<ctll::anything, atom_expr_impl_h<P>>;
        template<typename P> static constexpr auto rule(atom_expr_impl_h<P>, digit_chars) -> ctll::push<ctll::anything, make_impl_h, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_impl_h<P>, not_digit_chars) -> ctll::push<make_has_impl_h, atom_expr2<P>>;

        // cyclic: 'r'
        // acyclic: 'r0'
        // ring size: 'r' NUMBER
        // FIXME: r1 & r2  = error??
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'r'>) -> ctll::push<ctll::anything, atom_expr_ring_size<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size<P>, ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size<P>, ctll::range<'1', '9'>) -> ctll::push<ctll::anything, start_number, atom_expr_ring_size2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size<P>, not_digit_chars) -> ctll::push<make_cyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_ring_size2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_size2<P>, not_digit_chars) -> ctll::push<make_ring_size, atom_expr2<P>>;

        // cyclic: 'x'
        // acyclic: 'x0'
        // ring connectivity: 'x' NUMBER
        // FIXME: x1 = error?
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'x'>) -> ctll::push<ctll::anything, atom_expr_ring_connectivity<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity<P>, ctll::term<'0'>) -> ctll::push<ctll::anything, make_acyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity<P>, ctll::range<'1', '9'>) -> ctll::push<ctll::anything, start_number, atom_expr_ring_connectivity2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity<P>, not_digit_chars) -> ctll::push<make_cyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_ring_connectivity2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_connectivity2<P>, not_digit_chars) -> ctll::push<make_ring_connectivity, atom_expr2<P>>;

        // charge: '-' | '--' | '---' | ... | '-' DIGIT
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'-'>) -> ctll::push<ctll::anything, start_charge, atom_expr_neg_charge<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge<P>, digit_chars) -> ctll::push<ctll::anything, make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge<P>, ctll::term<'-'>) -> ctll::push<ctll::anything, decrement_charge, atom_expr_neg_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge<P>, ctll::neg_set<'-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge2<P>, ctll::term<'-'>) -> ctll::push<ctll::anything, decrement_charge, atom_expr_neg_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_neg_charge2<P>, ctll::neg_set<'-'>) -> ctll::push<make_charge, atom_expr2<P>>;

        // charge: '+' | '++' | '+++' | ... | '+' DIGIT
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'+'>) -> ctll::push<ctll::anything, start_charge, atom_expr_pos_charge<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge<P>, digit_chars) -> ctll::push<ctll::anything, make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge<P>, ctll::term<'+'>) -> ctll::push<ctll::anything, increment_charge, atom_expr_pos_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge<P>, ctll::neg_set<'+', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_charge, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge2<P>, ctll::term<'+'>) -> ctll::push<ctll::anything, increment_charge, atom_expr_pos_charge2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_pos_charge2<P>, ctll::neg_set<'+'>) -> ctll::push<make_charge, atom_expr2<P>>;

        // atom class: ':' NUMBER
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<':'>) -> ctll::push<ctll::anything, atom_expr_class<P>>;
        template<typename P> static constexpr auto rule(atom_expr_class<P>, digit_chars) -> ctll::push<ctll::anything, start_number, atom_expr_class2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_class<P>, not_digit_chars) -> ctll::reject; //ctll::push<make_class, set_and_high, atom_expr2>;
        template<typename P> static constexpr auto rule(atom_expr_class2<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_class2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_class2<P>, not_digit_chars) -> ctll::push<make_class, atom_expr2<P>>;

        // symbol: 'U' | 'V' | 'W'
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'U', 'V', 'W'>) -> ctll::push<ctll::anything, make_aliphatic, atom_expr2<P>>;
        // symbol: 'b' | 'c' | 'n' | 'o' | 'p'
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'b', 'c', 'n', 'o', 'p'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2<P>>;
        // any atom: '*'
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'*'>) -> ctll::push<ctll::anything, make_any_atom, atom_expr2<P>>;

        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'A'>) -> ctll::push<ctll::anything, push_char, atom_exprA<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'B'>) -> ctll::push<ctll::anything, push_char, atom_exprB<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'C'>) -> ctll::push<ctll::anything, push_char, atom_exprC<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'D'>) -> ctll::push<ctll::anything, push_char, atom_exprD<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'E'>) -> ctll::push<ctll::anything, push_char, atom_exprE<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'F'>) -> ctll::push<ctll::anything, push_char, atom_exprF<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'G'>) -> ctll::push<ctll::anything, push_char, atom_exprG<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'H'>) -> ctll::push<ctll::anything, push_char, atom_exprH<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'I'>) -> ctll::push<ctll::anything, push_char, atom_exprI<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'K'>) -> ctll::push<ctll::anything, push_char, atom_exprK<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'L'>) -> ctll::push<ctll::anything, push_char, atom_exprL<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'M'>) -> ctll::push<ctll::anything, push_char, atom_exprM<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'N'>) -> ctll::push<ctll::anything, push_char, atom_exprN<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'O'>) -> ctll::push<ctll::anything, push_char, atom_exprO<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'P'>) -> ctll::push<ctll::anything, push_char, atom_exprP<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'R'>) -> ctll::push<ctll::anything, push_char, atom_exprR<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'S'>) -> ctll::push<ctll::anything, push_char, atom_exprS<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'T'>) -> ctll::push<ctll::anything, push_char, atom_exprT<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'X'>) -> ctll::push<ctll::anything, push_char, atom_exprX<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'Y'>) -> ctll::push<ctll::anything, push_char, atom_exprY<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'Z'>) -> ctll::push<ctll::anything, push_char, atom_exprZ<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'a'>) -> ctll::push<ctll::anything, push_char, atom_expr_a<P>>;
        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::set<'s'>) -> ctll::push<ctll::anything, push_char, atom_expr_s<P>>;

        template<typename P> using Symbol1 = ctll::push<pop_char, make_aliphatic, atom_expr2<P>>;
        template<typename P> using Symbol2 = ctll::push<ctll::anything, make_aliphatic, atom_expr2<P>>;

        // symbol: 'Ac' | 'Ag' | 'Al' | 'Am' | 'Ar' | 'As' | 'At' | 'Au'
        // any aliphatic: 'A'
        template<typename P> static constexpr auto rule(atom_exprA<P>,     ctll::set<'c', 'g', 'l', 'm', 'r', 's', 't', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprA<P>, ctll::neg_set<'c', 'g', 'l', 'm', 'r', 's', 't', 'u'>) -> ctll::push<pop_char, make_any_aliphatic, atom_expr2<P>>;
        // symbol: 'B' | 'Ba' | 'Be' | 'Bh' | 'Bi' | 'Bk' | 'Br'
        template<typename P> static constexpr auto rule(atom_exprB<P>,     ctll::set<'a', 'e', 'h', 'i', 'k', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprB<P>, ctll::neg_set<'a', 'e', 'h', 'i', 'k', 'r'>) -> Symbol1<P>;
        // C Ca Cd Ce Cf Cl Cm Co Cr Cs Cu
        template<typename P> static constexpr auto rule(atom_exprC<P>,     ctll::set<'a', 'd', 'e', 'f', 'l', 'm', 'o', 'r', 's', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprC<P>, ctll::neg_set<'a', 'd', 'e', 'f', 'l', 'm', 'o', 'r', 's', 'u'>) -> Symbol1<P>;
        // symbol: 'Db' | 'Ds' | 'Dy'
        // degree: 'D' | 'D' NUMBER
        template<typename P> static constexpr auto rule(atom_exprD<P>,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_degree<P>>;
        template<typename P> static constexpr auto rule(atom_exprD<P>,     ctll::set<'b', 's', 'y'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprD<P>, ctll::neg_set<'b', 's', 'y', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_degree, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_degree<P>,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, push_number, atom_expr_degree<P>>;
        template<typename P> static constexpr auto rule(atom_expr_degree<P>, ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_degree, atom_expr2<P>>;
        // symbol: 'Er' | 'Es' | 'Eu'
        template<typename P> static constexpr auto rule(atom_exprE<P>,     ctll::set<'r', 's', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprE<P>, ctll::neg_set<'r', 's', 'u'>) -> ctll::reject;
        // symbol: 'F' | 'Fe' | 'Fm' | 'Fr'
        template<typename P> static constexpr auto rule(atom_exprF<P>,     ctll::set<'e', 'm', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprF<P>, ctll::neg_set<'e', 'm', 'r'>) -> Symbol1<P>;
        // symbol: 'Ga' | 'Gd' | 'Ge'
        template<typename P> static constexpr auto rule(atom_exprG<P>,     ctll::set<'a', 'd', 'e'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprG<P>, ctll::neg_set<'a', 'd', 'e'>) -> ctll::reject;
        // symbol: 'H' | 'He' | 'Hf' | 'Hg' | 'Ho' | 'Hs'
        // total hydrogens: 'H' | 'H' DIGIT
        template<typename P> static constexpr auto rule(atom_exprH<P>,     ctll::set<'e', 'f', 'g', 'o', 's'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprH<P>,   ctll::range<'0', '9'>) -> ctll::list<ctll::anything, pop_char, make_total_h, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_exprH<P>, ctll::neg_set<'e', 'f', 'g', 'o', 's', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::list<pop_char, make_total_h, atom_expr2<P>>;
        // symbol: 'I' | 'In' | 'Ir'
        template<typename P> static constexpr auto rule(atom_exprI<P>,     ctll::set<'n', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprI<P>, ctll::neg_set<'n', 'r'>) -> Symbol1<P>;
        // symbol: 'K' | 'Kr'
        template<typename P> static constexpr auto rule(atom_exprK<P>,    ctll::term<'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprK<P>, ctll::neg_set<'r'>) -> Symbol1<P>;
        // La Li Lr Lu
        template<typename P> static constexpr auto rule(atom_exprL<P>,     ctll::set<'a', 'i', 'r', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprL<P>, ctll::neg_set<'a', 'i', 'r', 'u'>) -> ctll::reject;
        // Md Mg Mn Mo Mt
        template<typename P> static constexpr auto rule(atom_exprM<P>,     ctll::set<'d', 'g', 'n', 'o', 't'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprM<P>, ctll::neg_set<'d', 'g', 'n', 'o', 't'>) -> ctll::reject;
        // N Na Nb Nd Ne Ni No Np
        template<typename P> static constexpr auto rule(atom_exprN<P>,     ctll::set<'a', 'b', 'd', 'e', 'i', 'o', 'p'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprN<P>, ctll::neg_set<'a', 'b', 'd', 'e', 'i', 'o', 'p'>) -> Symbol1<P>;
        // O Os
        template<typename P> static constexpr auto rule(atom_exprO<P>,    ctll::term<'s'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprO<P>, ctll::neg_set<'s'>) -> Symbol1<P>;
        // P Pa Pb Pd Pm Po Pr Pt Pu
        template<typename P> static constexpr auto rule(atom_exprP<P>,     ctll::set<'a', 'b', 'd', 'm', 'o', 'r', 't', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprP<P>, ctll::neg_set<'a', 'b', 'd', 'm', 'o', 'r', 't', 'u'>) -> Symbol1<P>;
        // symbol: 'Ra' | 'Rb' | 'Re' | 'Rf' | 'Rg' | 'Rh' | 'Rn' | 'Ru'
        // cyclic: 'R'
        // acyclic: 'R0'
        // ring count: 'R' NUMBER
        template<typename P> static constexpr auto rule(atom_exprR<P>,   ctll::term<'0'>) -> ctll::push<ctll::anything, pop_char, make_acyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_exprR<P>,   ctll::range<'1', '9'>) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_ring_count<P>>;
        template<typename P> static constexpr auto rule(atom_exprR<P>,     ctll::set<'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprR<P>, ctll::neg_set<'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_cyclic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_count<P>,   ctll::range<'0', '9'>) -> ctll::push<ctll::anything, push_number, atom_expr_ring_count<P>>;
        template<typename P> static constexpr auto rule(atom_expr_ring_count<P>, ctll::neg_set<'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<make_ring_count, atom_expr2<P>>;
        // S Sb Sc Se Sg Si Sm Sn Sr
        template<typename P> static constexpr auto rule(atom_exprS<P>,     ctll::set<'b', 'c', 'e', 'g', 'i', 'm', 'n', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprS<P>, ctll::neg_set<'b', 'c', 'e', 'g', 'i', 'm', 'n', 'r'>) -> Symbol1<P>;
        // Ta Tb Tc Te Th Ti Tl Tm
        template<typename P> static constexpr auto rule(atom_exprT<P>,     ctll::set<'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprT<P>, ctll::neg_set<'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm'>) -> ctll::reject;
        // symbol: 'Xe'
        // connectivity: 'X' | 'X' NUMBER
        template<typename P> static constexpr auto rule(atom_exprX<P>, digit_chars) -> ctll::push<ctll::anything, pop_char, start_number, atom_expr_connectivity<P>>;
        template<typename P> static constexpr auto rule(atom_exprX<P>,    ctll::term<'e'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprX<P>, ctll::neg_set<'e', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>) -> ctll::push<pop_char, make_connectivity, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_connectivity<P>, digit_chars) -> ctll::push<ctll::anything, push_number, atom_expr_connectivity<P>>;
        template<typename P> static constexpr auto rule(atom_expr_connectivity<P>, not_digit_chars) -> ctll::push<make_connectivity, atom_expr2<P>>;
        // symbol: 'Y' | 'Yb'
        template<typename P> static constexpr auto rule(atom_exprY<P>,    ctll::term<'b'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprY<P>, ctll::neg_set<'b'>) -> Symbol1<P>;
        // symbol: 'Zn' | 'Zr'
        template<typename P> static constexpr auto rule(atom_exprZ<P>,     ctll::set<'n', 'r'>) -> Symbol2<P>;
        template<typename P> static constexpr auto rule(atom_exprZ<P>, ctll::neg_set<'n', 'r'>) -> ctll::reject;

        // 'a' | 'as'
        template<typename P> static constexpr auto rule(atom_expr_a<P>,     ctll::set<'s'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_a<P>, ctll::neg_set<'s'>) -> ctll::push<pop_char, make_any_aromatic, atom_expr2<P>>;
        // 's' | 'se'
        template<typename P> static constexpr auto rule(atom_expr_s<P>,     ctll::set<'e'>) -> ctll::push<ctll::anything, make_aromatic, atom_expr2<P>>;
        template<typename P> static constexpr auto rule(atom_expr_s<P>, ctll::neg_set<'e'>) -> ctll::push<pop_char, make_aromatic, atom_expr2<P>>;

        //
        // Recursive SMARTS
        //

        template<typename P> static constexpr auto rule(atom_expr2<P>, ctll::term<'$'>) -> ctll::push<ctll::anything, atom_expr_recursive<P>>;
        template<typename P> static constexpr auto rule(atom_expr_recursive<P>, ctll::term<'('>) -> ctll::push<ctll::anything, push_recursive, atom<ParenthesisPushRecursive<P>>>;
        template<typename P> static constexpr auto rule(atom_expr_recursive<P>, ctll::neg_set<'('>) -> ctll::reject;

        //
        // Bond expressions
        //

        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::set<'-', '=', '#', '$', ':', '~', '@'>) -> ctll::push<ctll::anything, make_bond_primitive, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::set<'/', '\\'>) -> ctll::push<ctll::anything, push_char, bond_expr2<P>>;

        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_bond_not, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<'&'>) -> ctll::push<ctll::anything, make_bond_and_high, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<','>) -> ctll::push<ctll::anything, make_bond_or, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::term<';'>) -> ctll::push<ctll::anything, make_bond_and_low, bond_expr<P>>;

        template<typename P> static constexpr auto rule(bond_expr<P>, ctll::neg_set<'-', '=', '#', '$', ':', '~', '@', '/', '\\', '!', '&', ',', ';'>) -> ctll::push<chain<P>>; // FIXME: ring bonds?

        template<typename P> static constexpr auto rule(bond_expr2<P>, ctll::term<'?'>) -> ctll::push<ctll::anything, make_bond_primitive, pop_char, bond_expr<P>>;
        template<typename P> static constexpr auto rule(bond_expr2<P>, ctll::neg_set<'?'>) -> ctll::push<make_bond_primitive, pop_char, bond_expr<P>>;

        //
        // Chain expressions
        //

        template<typename P> static constexpr auto rule(chain<P>, ctll::epsilon) -> ctll::epsilon;
        template<typename P> static constexpr auto rule(chain<P>, ctll::set<'-', '=', '#', '$', ':', '~', '@', '/', '\\'>) -> ctll::push<bond_expr<P>>;
        //template<typename P> static constexpr auto rule(chain, ctll::set<'/', '\\'>) -> ctll::push<ctll::anything, make_bond_primitive, chain_up_down<P>>;
        //template<typename P> static constexpr auto rule(chain_up_down, ctll::neg_set<'?'>) -> ctll::push<chain<P>>;
        //template<typename P> static constexpr auto rule(chain_up_down, ctll::term<'?'>) -> ctll::push<ctll::anything, make_bond_primitive, set_and_high, reset_not, atom<P>>;

        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'!'>) -> ctll::push<ctll::anything, make_bond_not, bond_expr<P>>;
        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'.'>) -> ctll::push<ctll::anything, reset_prev, chain<P>>;

        //
        // Parenthesis
        //

        // push_prev
        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'('>) -> ctll::push<ctll::anything, push_prev, chain<ParenthesisPushBranch<P>>>;

        template<typename P> static constexpr auto rule(chain<P>, ctll::term<')'>)
        {
            if constexpr (P::branch > 0)
                // pop_prev
                return ctll::push<ctll::anything, pop_prev, chain<ParenthesisPopBranch<P>>>{};
            else if constexpr (!ctll::empty(P::recursive))
                // pop_recursive
                return ctll::push<ctll::anything, pop_recursive, atom_expr2<ParenthesisPopRecursive<P>>>{};
            else if constexpr (P::component)
                // pop_component
                return ctll::push<ctll::anything, pop_prev, chain<Parenthesis<>>>{};
            else
                return ctll::reject{};
        }

        // ring bond: DIGIT | '%' DIGIT DIGIT
        template<typename P> static constexpr auto rule(chain<P>, digit_chars) -> ctll::push<ctll::anything, handle_ring_bond, chain<P>>;
        template<typename P> static constexpr auto rule(chain<P>, ctll::term<'%'>) -> ctll::push<ctll::anything, ring_bond<P>>;
        template<typename P> static constexpr auto rule(ring_bond<P>, digit_chars) -> ctll::push<ctll::anything, start_number, ring_bond2<P>>;
        template<typename P> static constexpr auto rule(ring_bond2<P>, digit_chars) -> ctll::push<ctll::anything, handle_ring_bond, chain<P>>;

        // FIXME: add ops  etc
        using not_chain_chars = ctll::neg_set<'-', '=', '#', '$', ':', '~', '@', '/', '\\', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'>;
        template<typename P> static constexpr auto rule(chain<P>, not_chain_chars) -> ctll::push<atom<P>>;

        // chain ::= branched_atom | chain branched_atom | chain bond branched_atom | chain dot branched_atom

    };

} // namespace ctsmarts

#include <ctll/list.hpp>
#include <ctll/fixed_string.hpp>

namespace Kitimar::CTSmarts {

    template<auto N>
    std::string toString(ctll::fixed_string<N> str)
    {
        return std::string{str.begin(), str.end()};
    }

    template<int Index, typename ...Ts>
    constexpr auto get(ctll::list<Ts...>)
    {
        return std::get<Index>(std::tuple<Ts...>());
    }

    // resize list by adding T at front
    template<int Size, typename T, typename ...Ts>
    constexpr auto resize(ctll::list<Ts...> = {})
    {
        if constexpr (ctll::size(ctll::list<Ts...>()) == Size)
            return ctll::list<Ts...>();
        else
            return resize<Size, T>(ctll::list<T, Ts...>());
    }

    template<int Index, typename U, typename T, typename ...Ts, typename ...Us>
    constexpr auto replace(ctll::list<T, Ts...>, ctll::list<Us...> = {})
    {
        if constexpr (Index)
            return replace<Index - 1, U>(ctll::list<Ts...>(), ctll::list<Us..., T>());
        else
            return ctll::list<Us..., U, Ts...>();
    }

    template<typename Separator, typename ...Ts, typename ...Us, typename ...Vs>
    constexpr auto split(ctll::list<Ts...> list, ctll::list<Us...> parts = ctll::empty_list(), ctll::list<Vs...> part = ctll::empty_list())
    {
        if constexpr (ctll::empty(list)) {
            static_assert(!ctll::empty(part));
            return ctll::list<Us..., ctll::list<Vs...>>();
        } else {
            auto [head, tail] = ctll::pop_and_get_front(list);
            if constexpr (std::is_same_v<decltype(head), Separator>) {
                static_assert(!ctll::empty(part));
                auto parts2 = ctll::list<Us..., ctll::list<Vs...>>();
                return split<Separator>(tail, parts2, ctll::empty_list());
            } else {
                auto part2 = ctll::list<Vs..., decltype(head)>();
                return split<Separator>(tail, parts, part2);
            }
        }
    }

    template<typename ...Ts, typename F>
    constexpr auto transform(ctll::list<Ts...> list, F &&f)
    {
        if constexpr (ctll::empty(list)) {
            return ctll::empty_list{};
        } else {
            auto [head, tail] = ctll::pop_and_get_front(list);
            return ctll::push_front(f(head), transform(tail, std::forward<F>(f)));
        }
    }

    /*
    constexpr auto unique(ctll::empty_list) noexcept
    {
        return ctll::empty_list{};
    }

    template<typename T, typename ...Ts>
    constexpr auto unique(ctll::list<T, Ts...> list) noexcept
    {
        if constexpr (ctll::exists_in(T{}, ctll::list<Ts...>{}))
            return unique(ctll::list<Ts...>{});
        else
            return ctll::push_front(T{}, unique(ctll::list<Ts...>{}));
    }
    */

    constexpr auto unique(auto list) noexcept
    {
        if constexpr (ctll::empty(list))
            return ctll::empty_list{};
        else {
            constexpr auto head = ctll::front(list);
            constexpr auto tail = unique(ctll::pop_front(list));
            if constexpr (ctll::exists_in(head, tail))
                return tail;
            else
                return ctll::push_front(head, tail);
        }
    }

    template<typename ...Ts, typename ...Us>
    constexpr auto merge(ctll::list<Ts...> a, ctll::list<Us...> b) noexcept
    {
        return unique(ctll::concat(a, b));
    }

    /*
    constexpr auto merge(auto a, auto b)
    {
        if constexpr (!ctll::size(a))
            return b;
        else {
            constexpr auto head = ctll::front(a);
            if constexpr (ctll::exists_in(head, b))
                return b;
            else
                return merge(ctll::pop_front(a), ctll::push_front(head, b));
        }
    }
    */

    template<int ...N>
    constexpr auto toArray(ctll::list<Number<N>...>)
    {
        return std::array<int, sizeof...(N)>({ N... });
    }

} // namespace Kitimar::CTSmarts

#include <type_traits>
#include <cctype>

namespace Kitimar::CTSmarts {

    template<int N, int AtomIndex, typename BondExpr>
    struct RingBondHelper
    {
        static constexpr inline auto n = N;
        static constexpr inline auto atomIndex = AtomIndex;
        static constexpr inline auto bondExpr = BondExpr();
    };

    template<auto N, auto AtomIndex>
    struct ClassHelper
    {
        static constexpr inline auto n = N;
        static constexpr inline auto atomIndex = AtomIndex;
    };

    // Operation tags
    struct NotTag {};
    struct AndHighTag {};
    struct OrTag {};
    struct AndLowTag {};

    // Error tags

    struct NoErrorTag {};

    struct EmptyBracketAtomTag {}; // '[]'
    struct OpenBracketAtomTag {}; // '[C'
    struct MissingSymbolTag {}; // '[13]' '[+]'
    struct InvalidAtomExprTag {};
    struct InvalidChiralTag {};
    struct InvalidClassTag {}; // '[C:]'
    struct InvalidAtomPrimitiveTag {}; // '[Q]'
    struct InvalidBondPrimitiveTag {}; // 'C^C'
    struct InvalidChiralValenceTag {};
    struct InvalidChiralHydrogenCountTag {};
    struct HydrogenHydrogenCountTag {}; // '[HH1]'

    struct UnmatchedBranchOpeningTag {}; // 'CC(C'
    struct UnmatchedBranchClosingTag {}; // 'CC)C'
    struct InvalidRingBondTag {}; // 'C%' 'C%1'
    struct ConflicingRingBondTag {}; // 'C-1CC=1'
    struct LoopRingBondTag {}; // 'C11'
    struct DuplicateRingBondTag {}; // 'C12CC12'

    struct LeadingDotTag {}; // '.C'
    struct TrailingDotTag {}; // 'C.'
    struct LeadingBondTag {}; // '-C'
    struct TrailingBondTag {}; // 'C-'

    struct MissingLeftOperandTag {}; // '[&C]'
    struct MissingRightOperandTag {}; // '[C&]' '[!]'

    template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    struct SmartsParams
    {
        using NextIndex = T1;
        using PrevIndex = T2;
        using AtomExpr = T3;
        using BondExpr = T4;
        using RingBonds = T5;
        using Classes = T6;
        using Error = T7;

        static constexpr inline auto nextIndex = NextIndex{}; // Number
        static constexpr inline auto prevIndex = PrevIndex{}; // ctll::list<Number>
        static constexpr inline auto atomExpr = AtomExpr{}; // ctll::list<T>
        static constexpr inline auto bondExpr = BondExpr{};  // ctll::list<T>
        static constexpr inline auto ringBonds = RingBonds{}; // ctll::list<RingBondHelper>
        static constexpr inline auto classes = Classes{}; // ctll::list<ClassHelper>
        static constexpr inline auto error = Error{};

        template<typename T> static consteval auto setNextIndex(T) noexcept { return SmartsParams<T, T2, T3, T4, T5, T6, T7>{}; }
        template<typename T> static consteval auto setPrevIndex(T) noexcept { return SmartsParams<T1, T, T3, T4, T5, T6, T7>{}; }
        template<typename T> static consteval auto setAtomExpr(T) noexcept  { return SmartsParams<T1, T2, T, T4, T5, T6, T7>{}; }
        template<typename T> static consteval auto setBondExpr(T) noexcept  { return SmartsParams<T1, T2, T3, T, T5, T6, T7>{}; }
        template<typename T> static consteval auto setRingBonds(T) noexcept { return SmartsParams<T1, T2, T3, T4, T, T6, T7>{}; }
        template<typename T> static consteval auto setClasses(T) noexcept   { return SmartsParams<T1, T2, T3, T4, T5, T, T7>{}; }
        template<typename T> static consteval auto setError(T) noexcept     { return SmartsParams<T1, T2, T3, T4, T5, T6, T>{}; }

        consteval SmartsParams() noexcept = default;

        constexpr auto pushPrevIndex() const
        {
            return setPrevIndex(ctll::push_front(ctll::front(prevIndex), prevIndex));
        }
        constexpr auto popPrevIndex() const
        {
            static_assert(!ctll::empty(prevIndex));
            return setPrevIndex(ctll::pop_front(prevIndex));
        }

        template<typename PrevIndexT>
        constexpr auto nextAtomHelper(PrevIndexT) const
        {
            return SmartsParams<Number<nextIndex.value + 1>, PrevIndexT, AtomExpr, ctll::empty_list, RingBonds, Classes, Error>{};
        }

        constexpr auto nextAtom() const
        {
            if constexpr (ctll::empty(prevIndex))
                return nextAtomHelper(ctll::push_front(nextIndex, prevIndex));
            else
                return nextAtomHelper(ctll::pop_front_and_push_front(nextIndex, prevIndex));
        }

    };

    template <typename Atoms = ctll::empty_list, typename Bonds = ctll::empty_list,
              typename ParamsT = SmartsParams<Number<0>, ctll::empty_list, ctll::empty_list, ctll::empty_list, ctll::empty_list, ctll::empty_list, NoErrorTag>,
              typename ParentT = ctll::_nothing>
    struct SmartsContext
    {
        using Params = ParamsT;
        using Parent = ParentT;

        static constexpr inline auto atoms = Atoms{};
        static constexpr inline auto bonds = Bonds{};
        static constexpr inline auto params = ParamsT{};
        static constexpr inline auto parent = ParentT{};

        static constexpr inline auto valid = !ctll::size(params.ringBonds);

        constexpr SmartsContext() noexcept {}
        constexpr SmartsContext(Atoms, Bonds, Params, Parent) noexcept {}
    };

    template<auto V>
    constexpr auto bondPrimitive(ctll::term<V>)
    {
        if constexpr (V == '-') return BondOrder<1>();
        else if constexpr (V == '=') return BondOrder<2>();
        else if constexpr (V == '#') return BondOrder<3>();
        else if constexpr (V == '$') return BondOrder<4>();
        else if constexpr (V == ':') return AromaticBond();
        else if constexpr (V == '~') return AnyBond();
        else if constexpr (V == '@') return RingBond();
        else if constexpr (V == '/') return UpBond();
        else if constexpr (V == '\\') return DownBond();
        else return ImplicitBond();//static_assert(false);
    }

    template<auto V>
    constexpr auto termAtomicNumber(ctll::term<V>)
    {
        switch (V) {
            case 'H': return 1;
            case 'B': case 'b': return 5;
            case 'C': case 'c': return 6;
            case 'N': case 'n': return 7;
            case 'O': case 'o': return 8;
            case 'F': return 9;
            case 'P': case 'p': return 15;
            case 'S': case 's': return 16;
            case 'K': return 19;
            case 'V': return 23;
            case 'Y': return 39;
            case 'I': return 53;
            case 'W': return 74;
            case 'U': return 92;
            default: return -1;
        }
    }

    template<auto C, auto V>
    constexpr auto symbolAtomicNumber(Char<C>, ctll::term<V>)
    {
        switch (C) {
            case 'a':
                return V == 's' ? 33 : -1;
            case 's':
                return V == 'e' ? 34 : -1;
            case 'A':
                switch (V) {
                    case 'c': return 89;
                    case 'g': return 47;
                    case 'l': return 13;
                    case 'm': return 95;
                    case 'r': return 18;
                    case 's': return 33;
                    case 't': return 85;
                    case 'u': return 79;
                    default: break;
                }
                break;
            case 'B':
                switch (V) {
                    case 'a': return 56;
                    case 'e': return 4;
                    case 'i': return 83;
                    case 'k': return 97;
                    case 'r': return 35;
                    default: break;
                }
                break;
            case 'C':
                switch (V) {
                    case 'a': return 20;
                    case 'd': return 48;
                    case 'e': return 58;
                    case 'f': return 98;
                    case 'l': return 17;
                    case 'm': return 96;
                    case 'o': return 27;
                    case 'r': return 24;
                    case 's': return 55;
                    case 'u': return 29;
                    default: break;
                }
                break;
            case 'D':
                switch (V) {
                    case 'y': return 66;
                    default: break;
                }
                break;
            case 'E':
                switch (V) {
                    case 'r': return 68;
                    case 's': return 99;
                    case 'u': return 63;
                    default: break;
                }
                break;
            case 'F':
                switch (V) {
                    case 'e': return 26;
                    case 'm': return 100;
                    case 'r': return 87;
                    default: break;
                }
                break;
            case 'G':
                switch (V) {
                    case 'a': return 31;
                    case 'd': return 64;
                    case 'e': return 32;
                    default: break;
                }
                break;
            case 'H':
                switch (V) {
                    case 'e': return 2;
                    case 'f': return 72;
                    case 'g': return 80;
                    case 'o': return 67;
                    default: break;
                }
                break;
            case 'I':
                switch (V) {
                    case 'n': return 49;
                    case 'r': return 77;
                    default: break;
                }
                break;
            case 'K':
                switch (V) {
                    case 'r': return 36;
                    default: break;
                }
                break;
            case 'L':
                switch (V) {
                    case 'a': return 57;
                    case 'i': return 3;
                    case 'r': return 103;
                    case 'u': return 71;
                    default: break;
                }
                break;
            case 'M':
                switch (V) {
                    case 'd': return 101;
                    case 'g': return 12;
                    case 'n': return 25;
                    case 'o': return 42;
                    default: break;
                }
                break;
            case 'N':
                switch (V) {
                    case 'a': return 11;
                    case 'b': return 41;
                    case 'd': return 60;
                    case 'e': return 10;
                    case 'i': return 28;
                    case 'o': return 102;
                    case 'p': return 93;
                    default: break;
                }
                break;
            case 'O':
                switch (V) {
                    case 's': return 76;
                    default: break;
                }
                break;
            case 'P':
                switch (V) {
                    case 'a': return 91;
                    case 'b': return 82;
                    case 'd': return 46;
                    case 'm': return 61;
                    case 'o': return 84;
                    case 'r': return 59;
                    case 't': return 78;
                    case 'u': return 94;
                    default: break;
                }
                break;
            case 'R':
                switch (V) {
                    case 'a': return 88;
                    case 'b': return 37;
                    case 'e': return 75;
                    case 'h': return 45;
                    case 'n': return 86;
                    case 'u': return 44;
                    default: break;
                }
                break;
            case 'S':
                switch (V) {
                    case 'b': return 51;
                    case 'c': return 21;
                    case 'e': return 34;
                    case 'i': return 14;
                    case 'm': return 62;
                    case 'n': return 50;
                    case 'r': return 38;
                    default: break;
                }
                break;
            case 'T':
                switch (V) {
                    case 'a': return 73;
                    case 'b': return 65;
                    case 'c': return 43;
                    case 'e': return 52;
                    case 'h': return 90;
                    case 'i': return 22;
                    case 'l': return 81;
                    case 'm': return 69;
                    default: break;
                }
                break;
            case 'X':
                switch (V) {
                    case 'e': return 54;
                    default: break;
                }
                break;
            case 'Y':
                switch (V) {
                    case 'b': return 70;
                    default: break;
                }
                break;
            case 'Z':
                switch (V) {
                    case 'n': return 30;
                    case 'r': return 40;
                    default: break;
                }
                break;
        }
    }

    struct SmartsActions
    {

        //
        // Helper actions for parsing strings & numbers
        //

        // push_char
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_char, ctll::term<V>, Context ctx)
        {
            auto atoms = ctll::push_front(Char<V>(), ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        // pop_char
        template <auto V, auto C, typename ... Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::pop_char, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return SmartsContext{ctll::list<Ts...>(), ctx.bonds, ctx.params, ctx.parent};
        }

        // start_number
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::start_number, ctll::term<V>, Context ctx)
        {
            auto atoms = ctll::push_front(Number<V - '0'>(), ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        // push_number
        template <auto V, auto N, typename ... Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::push_number, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return SmartsContext{ctll::list<Number<10 * N + V - '0'>, Ts...>(), ctx.bonds, ctx.params, ctx.parent};
        }

        //
        // Operator helpers to convert infix expression list with primitives
        // and operator tags to the corresponding AST
        //

        template<typename Expr>
        static constexpr auto isOperation(Expr) { return false; }
        static constexpr auto isOperation(NotTag) { return true; }
        static constexpr auto isOperation(AndHighTag) { return true; }
        static constexpr auto isOperation(OrTag) { return true; }
        static constexpr auto isOperation(AndLowTag) { return true; }

        template<typename ...Expr>
        static constexpr auto isImplicitAnd(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::empty(expr))
                return false;
            else
                return !isOperation(ctll::front(expr));
        }

        static constexpr auto pushExpr(auto list, auto expr)
        {
            if constexpr (isImplicitAnd(list))
                return ctll::push_front(expr, ctll::push_front(AndHighTag(), list));
            else
                return ctll::push_front(expr, list);
        }

        template<typename ...Expr>
        static constexpr auto makeNot(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) % 2 == 0)
                return Not(ctll::front(ctll::rotate(expr)));
            else
                return ctll::front(ctll::rotate(expr));
        }

        template<typename ...Expr>
        static constexpr auto makeAndHighHelper(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) == 1)
                return makeNot(ctll::front(expr));
            else
                return And(makeNot(Expr())...);
        }

        static constexpr auto makeAndHigh(auto expr)
        {
            return makeAndHighHelper(split<AndHighTag>(expr));
        }

        template<typename ...Expr>
        static constexpr auto makeOrHelper(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) == 1)
                return makeAndHigh(ctll::front(expr));
            else
                return Or(makeAndHigh(Expr())...);
        }

        static constexpr auto makeOr(auto expr)
        {
            return makeOrHelper(split<OrTag>(expr));
        }

        template<typename ...Expr>
        static constexpr auto makeAndLowHelper(ctll::list<Expr...> expr)
        {
            if constexpr (ctll::size(expr) == 1)
                return makeOr(ctll::front(expr));
            else
                return And(makeOr(Expr())...);
        }

        static constexpr auto makeAndLow(auto expr)
        {
            return makeAndLowHelper(split<AndLowTag>(expr));
        }

        static constexpr auto makeAtomAST(auto expr)
        {
            constexpr auto expr2 = makeAndLow(expr);
            if constexpr (SmartsActions::isTotalHExpr(expr2))
                return SmartsActions::toTotalHExpr(expr2);
            else
                return expr2;
        }

        static constexpr auto makeBondAST(auto expr)
        {
            if constexpr (ctll::empty(expr))
                return ImplicitBond();
            else
                return makeAndLow(expr);
        }

        //
        // Transition to next atom:
        // - convert atomExpr to AST and add it to atoms
        // - create bond if needed
        //

        template<typename Atoms, typename Bonds, typename Params, typename Parent>
        static constexpr auto makeBond(Atoms atoms, Bonds bonds, Params params, Parent parent)
        {
            if constexpr (ctll::empty(params.prevIndex)) {
                return SmartsContext{atoms, bonds, params.nextAtom(), parent};
            } else {
                constexpr auto prevIndex = ctll::front(params.prevIndex).value;
                auto expr = makeBondAST(ctll::rotate(params.bondExpr));
                auto bond = Bond<ctll::size(bonds), prevIndex, params.nextIndex(), decltype(expr)>();
                auto bonds2 = ctll::push_front(bond, bonds);
                return SmartsContext{atoms, bonds2, params.nextAtom(), parent};
            }
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::next_atom, ctll::term<V> term, Context ctx)
        {
            auto expr = makeAtomAST(ctll::rotate(ctx.params.atomExpr));
            auto atom = Atom<ctll::size(ctx.atoms), decltype(expr)>{};
            auto atoms = ctll::push_front(atom, ctx.atoms);
            return makeBond(atoms, ctx.bonds, ctx.params.setAtomExpr(ctll::empty_list()), ctx.parent); // FIXME empty atomExpr in nextAtom
        }

        //
        // Atom primitives
        //

        template<typename Context, typename Atoms, typename Leaf>
        static constexpr auto pushAtomExpr(Context ctx, Atoms atoms, Leaf leaf)
        {
            auto atomExpr = pushExpr(ctx.params.atomExpr, leaf);
            auto params = ctx.params.setAtomExpr(atomExpr);
            return SmartsContext{atoms, ctx.bonds, params, ctx.parent};
        }

        // make_any_atom
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_atom, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAtom());
        }

        // make_any_aliphatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_aliphatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAliphatic());
        }

        // make_any_aromatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_any_aromatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AnyAromatic());
        }

        // make_aliphatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_aliphatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AliphaticAtom<termAtomicNumber(term)>());
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_aliphatic, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            auto expr = AliphaticAtom<symbolAtomicNumber(Char<C>(), term)>();
            return pushAtomExpr(ctx, ctll::list<Ts...>(), expr);
        }

        // make_aromatic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_aromatic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, AromaticAtom<termAtomicNumber(term)>());
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_aromatic, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            auto expr = AromaticAtom<symbolAtomicNumber(Char<C>(), term)>();
            return pushAtomExpr(ctx, ctll::list<Ts...>(), expr);
        }

        // make_isotope
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_isotope, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Isotope<N>());
        }

        // make_element
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_element, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Element<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_element, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Element<N>());
        }

        // make_degree
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_degree, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Degree<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_degree, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Degree<N>());
        }

        // make_valence
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_valence, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Valence<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_valence, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Valence<N>());
        }

        // make_connectivity
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_connectivity, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Connectivity<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_connectivity, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), Connectivity<N>());
        }

        // make_cyclic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_cyclic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Cyclic());
        }

        // make_acyclic
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_acyclic, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, Acyclic());
        }

        // make_ring_count
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_count, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingCount<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_ring_count, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingCount<N>());
        }

        // make_ring_size
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_size, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingSize<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_ring_size, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingSize<N>());
        }

        // make_ring_connectivity
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_ring_connectivity, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, RingConnectivity<1>());
        }
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_ring_connectivity, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return pushAtomExpr(ctx, ctll::list<Ts...>(), RingConnectivity<N>());
        }

        // make_class
        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_class, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            auto cls = ctll::push_front(ClassHelper<N, ctx.params.nextIndex()>(), ctx.params.classes);
            return SmartsContext{ctll::list<Ts...>(), ctx.bonds, ctx.params.setClasses(cls), ctx.parent};
        }

        // make_has_impl_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_has_impl_h, ctll::term<V> term, Context ctx)
        {
            return pushAtomExpr(ctx, ctx.atoms, HasImplicitH{});
        }

        // make_impl_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_impl_h, ctll::term<V> term, Context ctx)
        {
            constexpr auto value = V == 'h' ? 1 : V - '0';
            return pushAtomExpr(ctx, ctx.atoms, ImplicitH<value>());
        }

        // aliphatic, aromatic, '*', 'A' or 'a'
        template<typename Expr>
        static constexpr bool isTotalHExpr(Expr) { return false; }
        static constexpr bool isTotalHExpr(AnyAtom) { return true; }
        static constexpr bool isTotalHExpr(AnyAliphatic) { return true; }
        static constexpr bool isTotalHExpr(AnyAromatic) { return true; }
        template<int AtomicNumber>
        static constexpr bool isTotalHExpr(AliphaticAtom<AtomicNumber>) { return true; }
        static constexpr bool isTotalHExpr(AliphaticAtom<1>) { return false; }
        template<int AtomicNumber>
        static constexpr bool isTotalHExpr(AromaticAtom<AtomicNumber>) { return true; }
        template<int AtomicNumber>
        static constexpr bool isTotalHExpr(Element<AtomicNumber>) { return true; }
        // FIXME: operations
        template<typename Expr>
        static constexpr bool isTotalHExpr(Not<Expr>) { return isTotalHExpr(Expr()); }
        template<typename ...Expr>
        static constexpr bool isTotalHExpr(And<Expr...>) { return (isTotalHExpr(Expr()) || ...); }
        template<typename ...Expr>
        static constexpr bool isTotalHExpr(Or<Expr...>) { return (isTotalHExpr(Expr()) || ...); }

        static constexpr bool isTotalH(auto atomExpr)
        {
            if constexpr (ctll::empty(atomExpr)) {
                return false;
            } else {
                auto [head, tail] = ctll::pop_and_get_front(atomExpr);
                if (isTotalHExpr(head))
                    return true;
                return isTotalH(tail);
            }
        }

        template<typename Expr>
        static constexpr auto toTotalHExpr(Expr) { return Expr{}; }
        static constexpr auto toTotalHExpr(AliphaticAtom<1>) { return TotalH<1>{}; }
        template<typename Expr>
        static constexpr auto toTotalHExpr(Not<Expr>) { return Not<decltype(toTotalHExpr(Expr{}))>{}; }

        template<typename ...Expr>
        static constexpr auto toTotalHExpr(And<Expr...> op) { return And(toTotalHExpr(op.expr)); }
        template<typename ...Expr>
        static constexpr auto toTotalHExpr(Or<Expr...> op) { return Or(toTotalHExpr(op.expr)); }

        template<typename ...Expr>
        static constexpr auto toTotalHExpr(ctll::list<Expr...> expr)
        {
            return transform(expr, [] <typename Expr2> (Expr2 expr2) {
                if constexpr (std::is_same_v<Expr2, AliphaticAtom<1>>)
                    return TotalH<1>{};
                else
                    return toTotalHExpr(expr2);
            });
        }

        // make_total_h
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_total_h, ctll::term<V> term, Context ctx)
        {
            if constexpr (V == 'H')
                // Will be replaced by TotalH later if needed
                return pushAtomExpr(ctx, ctx.atoms, AliphaticAtom<1>());
            else
                return pushAtomExpr(ctx, ctx.atoms, TotalH<V - '0'>());
        }

        // start_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::start_charge, ctll::term<V> term, Context ctx)
        {
            constexpr auto value = V == '+' ? 1 : -1;
            auto atoms = ctll::push_front(Charge<value>{}, ctx.atoms);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        template <int Value, typename Context>
        static constexpr auto add_charge(Context ctx)
        {
            auto [charge, tail] = ctll::pop_and_get_front(ctx.atoms);
            auto atoms = ctll::push_front(Charge<charge.value + Value>{}, tail);
            return SmartsContext{atoms, ctx.bonds, ctx.params, ctx.parent};
        }

        // increment_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::increment_charge, ctll::term<V> term, Context ctx)
        {
            return add_charge<1>(ctx);
        }

        // decrement_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::decrement_charge, ctll::term<V> term, Context ctx)
        {
            return add_charge<-1>(ctx);
        }

        // make_charge
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_charge, ctll::term<V> term, Context ctx)
        {
            auto [charge, atoms] = ctll::pop_and_get_front(ctx.atoms);
            constexpr bool isDigit = V >= '0' && V <= '9';
            if constexpr (isDigit)
                return pushAtomExpr(ctx, atoms, Charge<charge.value * (V - '0')>{});
            else
                return pushAtomExpr(ctx, atoms, charge);
        }

        //
        // Operations
        //

        // make_atom_not
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_not, ctll::term<V> term, Context ctx)
        {
            auto atomExpr = pushExpr(ctx.params.atomExpr, NotTag());
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setAtomExpr(atomExpr), ctx.parent};
        }

        static constexpr auto makeAtomOp(auto ctx, auto op)
        {
            auto atomExpr = ctll::push_front(op, ctx.params.atomExpr);
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setAtomExpr(atomExpr), ctx.parent};
        }

        // make_atom_and_high
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_and_high, ctll::term<V> term, Context ctx)
        {
            return makeAtomOp(ctx, AndHighTag());
        }

        // make_atom_or
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_or, ctll::term<V> term, Context ctx)
        {
            return makeAtomOp(ctx, OrTag());
        }

        // set_atom_and_low
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_atom_and_low, ctll::term<V> term, Context ctx)
        {
            return makeAtomOp(ctx, AndLowTag());
        }

        // make_bond_not
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_not, ctll::term<V> term, Context ctx)
        {
            auto bondExpr = pushExpr(ctx.params.bondExpr, NotTag());
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(bondExpr), ctx.parent};
        }

        static constexpr auto makeBondOp(auto ctx, auto op)
        {
            auto bondExpr = ctll::push_front(op, ctx.params.bondExpr);
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(bondExpr), ctx.parent};
        }

        // make_bond_and_high
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_and_high, ctll::term<V> term, Context ctx)
        {
            return makeBondOp(ctx, AndHighTag());
        }

        // make_bond_or
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_or, ctll::term<V> term, Context ctx)
        {
            return makeBondOp(ctx, OrTag());
        }

        // set_bond_and_low
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_and_low, ctll::term<V> term, Context ctx)
        {
            return makeBondOp(ctx, AndLowTag());
        }

        //
        // Branches
        //

        // push_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_prev, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.pushPrevIndex(), ctx.parent};
        }

        // pop_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::pop_prev, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.popPrevIndex(), ctx.parent};
        }

        // reset_prev
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::reset_prev, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.popPrevIndex(), ctx.parent};
        }

        //
        // Recursive
        //

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::push_recursive, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{{}, {}, {}, ctx};
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::pop_recursive, ctll::term<V> term, Context ctx)
        {
            auto expr = BasicSmarts{ctll::rotate(Context::atoms), ctll::rotate(Context::bonds), Context::params.classes};
            return pushAtomExpr(ctx.parent, ctx.parent.atoms, expr);
        }

        //
        // Bond primitives
        //

        template<typename Context, typename Atoms, typename Leaf>
        static constexpr auto pushBondExpr(Context ctx, Atoms atoms, Leaf leaf)
        {
            auto bondExpr = pushExpr(ctx.params.bondExpr, leaf);
            auto params = ctx.params.setBondExpr(bondExpr);
            return SmartsContext{atoms, ctx.bonds, params, ctx.parent};
        }

        // bond_primitive
        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::make_bond_primitive, ctll::term<V> term, Context ctx)
        {
            static_assert(!ctll::empty(ctx.params.prevIndex));
            return pushBondExpr(ctx, ctx.atoms, bondPrimitive(term));
        }
        template <auto V, auto C, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::make_bond_primitive, ctll::term<V> term, SmartsContext<ctll::list<Char<C>, Ts...>, Bonds, Params, Parent> ctx)
        {
            static_assert(!ctll::empty(ctx.params.prevIndex));
            if constexpr (V == '?') {
                auto expr = leafBondExpr<ctx.params.operation, ctx.params.notSet>(ctx.params.bondExpr, UpOrDownBond());
                return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr), ctx.parent};
            } else {
                auto expr = leafBondExpr<ctx.params.operation, ctx.params.notSet>(ctx.params.bondExpr, bondPrimitive(ctll::term<C>()));
                return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setBondExpr(expr), ctx.parent};
            }
        }

        //
        // Ring bonds
        //

        template<int N, typename ...RingBonds>
        static constexpr auto findRingBond(ctll::list<RingBonds...> ringBonds)
        {
            if constexpr (ctll::empty(ringBonds)) {
                return ctll::_nothing();
            } else {
                auto rb = ctll::front(ringBonds);
                if constexpr (std::is_same_v<decltype(rb), ctll::_nothing>)
                    return findRingBond<N>(ctll::pop_front(ringBonds));
                if constexpr (rb.n == N)
                    return rb;
                else
                    return findRingBond<N>(ctll::pop_front(ringBonds));
            }
        }

        template<int Index, int Source, int Target, typename Expr1, typename Expr2>
        static constexpr auto makeRingBond(Expr1, Expr2)
        {
            if constexpr (std::is_same_v<Expr1, Expr2>)
                return std::make_tuple(Bond<Index, Source, Target, Expr1>{}, NoErrorTag{});
            else if constexpr (std::is_same_v<Expr1, ImplicitBond>)
                return std::make_tuple(Bond<Index, Source, Target, Expr2>{}, NoErrorTag{});
            else if constexpr (std::is_same_v<Expr2, ImplicitBond>)
                return std::make_tuple(Bond<Index, Source, Target, Expr1>{}, NoErrorTag{});
            else
                return std::make_tuple(Bond<Index, Source, Target, Expr1>{}, ConflicingRingBondTag{});
        }

        // ring_bond
        template <auto N, typename Context, typename Atoms>
        static constexpr auto handleRingBond(Context ctx, Atoms atoms)
        {
            constexpr auto atomIndex = ctll::front(ctx.params.prevIndex).value;
            auto rb = findRingBond<N>(ctx.params.ringBonds);
            if constexpr (std::is_same_v<decltype(rb), ctll::_nothing>) {
                auto rb2 = RingBondHelper<N, atomIndex, decltype(ctx.params.bondExpr)>();
                auto ringBonds = ctll::push_front(rb2, ctx.params.ringBonds);
                return SmartsContext{atoms, ctx.bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()), ctx.parent};
            } else {
                constexpr auto prevIndex = rb.atomIndex;
                constexpr auto ringBond = makeRingBond<ctll::size(ctx.bonds), atomIndex, prevIndex>(makeBondAST(ctll::rotate(rb.bondExpr)), makeBondAST(ctll::rotate(ctx.params.bondExpr)));
                constexpr auto bond = std::get<0>(ringBond);
                constexpr auto error = std::get<1>(ringBond);
                auto bonds = ctll::push_front(bond, ctx.bonds);
                auto ringBonds = ctll::remove_item(rb, ctx.params.ringBonds);
                if constexpr (std::is_same_v<NoErrorTag, decltype(error)>)
                    return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()), ctx.parent};
                else
                    return SmartsContext{atoms, bonds, ctx.params.setRingBonds(ringBonds).setBondExpr(ctll::empty_list()).setError(error), ctx.parent};
            }
        }

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::handle_ring_bond, ctll::term<V> term, Context ctx)
        {
            return handleRingBond<V - '0'>(ctx, ctx.atoms);
        }

        template <auto V, auto N, typename ...Ts, typename Bonds, typename Params, typename Parent>
        static constexpr auto apply(SmartsGrammar::handle_ring_bond, ctll::term<V> term, SmartsContext<ctll::list<Number<N>, Ts...>, Bonds, Params, Parent> ctx)
        {
            return handleRingBond<10 * N + V - '0'>(ctx, ctll::list<Ts...>());
        }

        //
        // Errors
        //

        template <auto V, typename Context>
        static constexpr auto apply(SmartsGrammar::error_empty_bracket, ctll::term<V> term, Context ctx)
        {
            return SmartsContext{ctx.atoms, ctx.bonds, ctx.params.setError(EmptyBracketAtomTag()), ctx.parent};
        }

    };

} // namespace ctsmarts

#include <ctll/parser.hpp>

#include <ranges>
#include <algorithm>

namespace Kitimar::CTSmarts {

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
        constexpr auto mapping = captureMappingHelper<1>(smarts.classes, ctll::empty_list());
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
            switch (vertexDegree.data[i]) {
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

    template <ctll::fixed_string SMARTS, bool IgnoreInvalid = false,
              typename Result = ctll::parser<SmartsGrammar, SMARTS, SmartsActions>::template output<SmartsContext<>>,
              typename Context = Result::output_type>
    struct Smarts : BasicSmarts<decltype(ctll::rotate(Context::atoms)), decltype(ctll::rotate(Context::bonds)), std::remove_const_t<decltype(Context::params.classes)>>
    {
        static constexpr inline auto smarts = SMARTS;
        static constexpr inline auto context = Context{};
        static constexpr inline auto valid = Result::is_correct;
        static constexpr inline auto position = Result::position;

        static constexpr auto input()
        {
            auto str = SMARTS | std::views::transform([] (auto c) { return static_cast<char>(c); });
            return std::string(str.begin(), str.end());
        }

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

        static_assert(IgnoreInvalid || valid);
        static_assert(IgnoreInvalid || std::is_same_v<const NoError, decltype(error)>); // FIXME
        //static_assert(IgnoreInvalid || error); // FIXME

    };

} // namespace ctsmarts

#include <array>
#include <vector>
#include <ranges>
#include <cassert>

#ifdef KITIMAR_WITH_IOSTREAM
#include <iostream>

template<std::integral I, auto N>
std::ostream& operator<<(std::ostream &os, const std::array<I, N> &map)
{
    os << "[";
    for (auto i = 0UL; i < map.size(); ++i)
        os << " " << map[i];
    os << " ]";
    return os;
}

template<std::integral I>
std::ostream& operator<<(std::ostream &os, const std::vector<I> &v)
{
    os << "[ ";
    for (auto i : v)
        os << i << " ";
    os << "]";
    return os;
}

#endif // KITIMAR_WITH_IOSTREAM

namespace Kitimar::CTSmarts {

    template<std::integral Index, int N>
    using IsomorphismMap = std::array<Index, N>;

    template<std::integral Index, int N>
    using IsomorphismMaps = std::vector<IsomorphismMap<Index, N>>;

    template<std::integral Index, int NumQueryAtoms>
    class LookupMap
    {
        public:
            static constexpr inline auto invalidIndex = static_cast<Index>(-1);
            using Map = IsomorphismMap<Index, NumQueryAtoms>;

            constexpr LookupMap() noexcept
            {
                m_map.fill(invalidIndex);
            }

            constexpr const Map& map() const noexcept
            {
                return m_map;
            }

            constexpr Index operator()(int queryAtomIndex) const noexcept
            {
                return m_map[queryAtomIndex];
            }

            constexpr bool contains(int queryAtomIndex, auto atomIndex) const noexcept
            {
                return m_map[queryAtomIndex] == atomIndex;
            }

            constexpr bool containsAtom(auto atomIndex) const noexcept
            {
                return std::ranges::find(m_map, atomIndex) != m_map.end();
            }

            constexpr bool containsQueryAtom(int queryAtomIndex) const noexcept
            {
                return m_map[queryAtomIndex] != invalidIndex;
            }

            constexpr void reset(auto numAtoms) const noexcept {}

            constexpr void add(int queryAtomIndex, auto atomIndex) noexcept
            {
                m_map[queryAtomIndex] = atomIndex;
            }

            constexpr void remove(int queryAtomIndex, auto atomIndex) noexcept
            {
                m_map[queryAtomIndex] = invalidIndex;
            }

        protected:
            Map m_map;
    };

    template<std::integral Index, int NumAtoms>
    class InverseMap : public LookupMap<Index, NumAtoms>
    {
        public:
            constexpr void reset(auto numAtoms)
            {
                m_mapped.clear();
                m_mapped.resize(numAtoms);
            }

            constexpr bool containsAtom(auto atomIndex) const noexcept
            {
                assert(atomIndex < m_mapped.size());
                return m_mapped[atomIndex];
            }

            constexpr void add(int queryAtomIndex, auto atomIndex) noexcept
            {
                LookupMap<Index, NumAtoms>::add(queryAtomIndex, atomIndex);
                assert(atomIndex < m_mapped.size());
                assert(!m_mapped[atomIndex]);
                m_mapped[atomIndex] = true;
            }

            constexpr void remove(int queryAtomIndex, auto atomIndex) noexcept
            {
                LookupMap<Index, NumAtoms>::remove(queryAtomIndex, atomIndex);
                assert(atomIndex < m_mapped.size());
                assert(m_mapped[atomIndex]);
                m_mapped[atomIndex] = false;
            }

        private:
            std::vector<uint8_t> m_mapped;
    };

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    template<typename VertexDegree, typename DfsBondList>
    struct SmartsQuery
    {
        static constexpr inline auto degrees = VertexDegree{}.data;
        static constexpr inline auto bonds = DfsBondList{}.data;

        consteval SmartsQuery() noexcept = default;
        consteval SmartsQuery(VertexDegree, DfsBondList) noexcept {};
    };

    enum class SeedType {
        Molecule, // start from any atom/bond
        Atom, // start from first atom
        Bond // start from first bond
    };

    template<SeedType T>
    using SeedTypeTag = std::integral_constant<SeedType, T>;

    struct NoOptimizer
    {
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = EdgeList{smarts};
            constexpr auto vertexDegree = VertexDegree{smarts, edgeList};
            constexpr auto incidentList = IncidentList{smarts, edgeList, vertexDegree};
            constexpr auto dfsEdges = DfsEdgeList{smarts, incidentList};
            constexpr auto cycleMembership = CycleMembership{smarts, incidentList};
            constexpr auto dfsBonds = DfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

    // FIXME: optimize atom expressions...
    struct FullOptimizer
    {
        template<typename SmartsT, SeedType SeedT>
        static constexpr auto create(SmartsT, SeedTypeTag<SeedT>) noexcept
        {
            constexpr auto smarts = SmartsT{};
            constexpr auto edgeList = EdgeList{smarts};
            constexpr auto vertexDegree = VertexDegree{smarts, edgeList};
            constexpr auto incidentList = IncidentList{smarts, edgeList, vertexDegree};
            constexpr auto atomFreq = AtomFrequency{smarts};
            constexpr auto optimizedIncidentList = OptimizeIncidentList{smarts, incidentList, atomFreq, std::false_type{}};
            constexpr auto sourceIndex = (SeedT == SeedType::Molecule) ? std::ranges::min_element(atomFreq.data) - std::begin(atomFreq.data) : 0;
            constexpr auto dfsEdges = DfsEdgeList{smarts, optimizedIncidentList, Number<sourceIndex>{}};
            constexpr auto cycleMembership = CycleMembership{smarts, optimizedIncidentList};
            constexpr auto dfsBonds = DfsBondList{smarts, dfsEdges, cycleMembership};
            return SmartsQuery{vertexDegree, dfsBonds};
        }
    };

}

namespace Kitimar::CTSmarts {

    enum class DefaultImplicitH
    {
        AtLeastOne, // Daylight, OpenSMARTS, RDKit
        ExactlyOne // OpenBabel
    };

    enum class Specialize : int
    {
        None,
        Atom  = 1,
        Bond  = 2,
        Chain = 4,
        Star  = 8,
        All = Atom | Bond | Chain | Star
    };

    inline constexpr bool operator&(Specialize lhs, Specialize rhs) noexcept
    {
        return static_cast<int>(lhs) & static_cast<int>(rhs);
    }

    template<DefaultImplicitH H, Specialize S, typename OptimizerT, template<std::integral, int N> class MapT>
    struct Config
    {
        // Matching behavior
        static constexpr inline auto defaultImplicitH = H;

        // Performance
        static constexpr inline auto specialize = S;

        using Optimizer = OptimizerT;

        template<std::integral Index, int N>
        using Map = MapT<Index, N>;

        static constexpr auto transformSmarts(auto smarts) noexcept
        {
            if constexpr (H == DefaultImplicitH::AtLeastOne)
                return smarts;
            else
                return BasicSmarts<decltype(replaceExpr(HasImplicitH{}, ImplicitH<1>{}, smarts.atoms)), decltype(smarts.bonds), decltype(smarts.classes)>{};
        }
    };

    using DefaultConfig = Config<DefaultImplicitH::AtLeastOne, Specialize::All, FullOptimizer, InverseMap>;

    using NoOptimizeConfig = Config<DefaultImplicitH::AtLeastOne, Specialize::None, NoOptimizer, InverseMap>;

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Operators
    //

    template<typename Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Not<Expr>)
    {
        return !matchAtomExpr(mol, atom, Expr{});
    }

    template<typename ...Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Or<Expr...>)
    {
        return (matchAtomExpr(mol, atom, Expr{}) || ...);
    }

    template<typename ...Expr>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, And<Expr...>)
    {
        return (matchAtomExpr(mol, atom, Expr{}) && ...);
    }

    template<typename Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, Not<Expr>)
    {
        return !matchBondExpr(mol, bond, Expr{});
    }

    template<typename ...Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, Or<Expr...>)
    {
        return (matchBondExpr(mol, bond, Expr{}) || ...);
    }

    template<typename ...Expr>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, And<Expr...>)
    {
        return (matchBondExpr(mol, bond, Expr{}) && ...);
    }

    //
    // Atoms
    //

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AnyAtom)
    {
        return true;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AnyAromatic)
    {
        return is_aromatic_atom(mol, atom);
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AnyAliphatic)
    {
        return !is_aromatic_atom(mol, atom);
    }

    template<int Element>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AliphaticAtom<Element>)
    {
        return get_element(mol, atom) == Element && !is_aromatic_atom(mol, atom);
    }

    template<int Element>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, AromaticAtom<Element>)
    {
        return get_element(mol, atom) == Element && is_aromatic_atom(mol, atom);
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Isotope<N>)
    {
        return get_isotope(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Element<N>)
    {
        return get_element(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Degree<N>)
    {
        return get_degree(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Valence<N>)
    {
        return get_valence(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Connectivity<N>)
    {
        return get_degree(mol, atom) + get_implicit_hydrogens(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, TotalH<N>)
    {
        return get_total_hydrogens(mol, atom) == N;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, HasImplicitH)
    {
        return get_implicit_hydrogens(mol, atom) > 0;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, ImplicitH<N>)
    {
        return get_implicit_hydrogens(mol, atom) == N;
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Cyclic)
    {
        return is_ring_atom(mol, atom);
    }

    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Acyclic)
    {
        return !is_ring_atom(mol, atom);
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, RingCount<N>)
    {
        return get_ring_count(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, RingSize<N>)
    {
        return is_in_ring_size(mol, atom, N);
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, RingConnectivity<N>)
    {
        return get_ring_degree(mol, atom) == N;
    }

    template<int N>
    constexpr bool matchAtomExpr(const auto &mol, const auto &atom, Charge<N>)
    {
        return get_charge(mol, atom) == N;
    }

    //
    // Recursive
    //

    namespace impl {
        template<typename SmartsT, typename Config, Molecule::Molecule Mol>
        constexpr bool match_atom(const Mol &mol, const auto &atom);
    }

    template<Molecule::Molecule Mol, typename Atoms, typename Bonds, typename Classes>
    constexpr bool matchAtomExpr(const Mol &mol, const auto &atom, BasicSmarts<Atoms, Bonds, Classes>)
    {
        return impl::match_atom<BasicSmarts<Atoms, Bonds, Classes>, DefaultConfig, Mol>(mol, atom);
    }

    //
    // Bonds
    //

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, ImplicitBond)
    {
        return get_order(mol, bond) == 1 || is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, AnyBond)
    {
        return true;
    }

    template<int Order>
    constexpr bool matchBondExpr(const auto &mol, const auto &bond, BondOrder<Order>)
    {
        return get_order(mol, bond) == Order && !is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, AromaticBond)
    {
        return is_aromatic_bond(mol, bond);
    }

    constexpr bool matchBondExpr(const auto &mol, const auto &bond, RingBond)
    {
        return is_ring_bond(mol, bond);
    }

} // namespace ctsmarts

namespace Kitimar::CTSmarts {

    struct UnconditionalFilter
    {
        static consteval bool enable(auto smarts) noexcept
        {
            return true;
        }
    };

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    struct NumAtomBondFilter : UnconditionalFilter
    {
        static constexpr bool reject(auto smarts, Molecule::Molecule auto &mol) noexcept
        {
            return num_atoms(mol) < smarts.numAtoms || num_bonds(mol) < smarts.numBonds;
        }
    };

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    template<typename Primitive, typename Expr>
    consteval auto expressionRequiresPrimitive(Expr, Primitive) noexcept
    {
        return std::is_same_v<Primitive, Expr>;
    }

    template<typename ...Expr>
    consteval auto expressionRequiresPrimitive(And<Expr...>, auto primitive) noexcept
    {
        return (expressionRequiresPrimitive(Expr{}, primitive) || ...);
    }

    template<typename ...Atoms>
    consteval auto smartsRequiresAtomPrimitiveHelper(ctll::list<Atoms...>, auto primitive) noexcept
    {
        return (expressionRequiresPrimitive(Atoms::expr, primitive) || ...);
    }

    consteval auto smartsRequiresAtomPrimitive(auto smarts, auto primitive) noexcept
    {
        return smartsRequiresAtomPrimitiveHelper(smarts.atoms, primitive);
    }

    template<typename MaxFrequency, int AtomicNumber>
    consteval auto smartsRequiresElement(auto smarts) noexcept
    {

        if constexpr (expressionFrequency(Element<AtomicNumber>{}) > MaxFrequency::value)
            return false;
        else
            return smartsRequiresAtomPrimitive(smarts, Element<AtomicNumber>{}) ||
                    smartsRequiresAtomPrimitive(smarts, AliphaticAtom<AtomicNumber>{}) ||
                    smartsRequiresAtomPrimitive(smarts, AromaticAtom<AtomicNumber>{});
    }

    template<typename Expr>
    consteval auto requiredExpressionPrimitives(Expr) noexcept
    {
        return ctll::list<Expr>{};
    }

    template<typename ...Expr>
    consteval auto requiredExpressionPrimitives(ctll::list<Expr...> expr) noexcept
    {
        if constexpr (ctll::empty(expr))
            return ctll::empty_list{};
        else {
            auto [head, tail] = ctll::pop_and_get_front(expr);
            auto tailPrimitives = requiredPrimitives(tail);
            if constexpr (ctll::exists_in(head, tail))
                return tailPrimitives;
            else
                return ctll::push_front(head, tailPrimitives);
        }
    }

    template<typename ...Expr>
    consteval auto requiredExpressionPrimitives(And<Expr...> op) noexcept
    {
        return requiredPrimitives(op.expr);
    }

    template<typename ...Expr>
    consteval auto requiredExpressionPrimitives(Or<Expr...>) noexcept
    {
        return ctll::empty_list{};
    }

    template<typename Expr>
    consteval auto requiredExpressionPrimitives(Not<Expr>) noexcept
    {
        return ctll::empty_list{};
    }

    template<typename ...Atoms>
    consteval auto requiredAtomPrimitives(ctll::list<Atoms...> atoms) noexcept
    {
        if constexpr (ctll::empty(atoms))
            return ctll::empty_list{};
        else
            return merge(requiredExpressionPrimitives(ctll::front(atoms).expr), requiredAtomPrimitives(ctll::pop_front(atoms)));
    }

    /*
    template<typename Expr>
    consteval auto requiredAtomPrimitivesHelper(ctll::list<Atoms...>) noexcept
    {

    }

    template<typename ...Atoms>
    consteval auto requiredAtomPrimitives(ctll::list<Atoms...>) noexcept
    {

    }
    */

    namespace impl {

        template<int N, int ...Values>
        consteval auto toArrayHelper(ctll::list<Number<Values>...> numbers) noexcept
        {
            if constexpr (ctll::empty(numbers))
                return std::array<int, N>{};
            else {
                auto [number, tail] = ctll::pop_and_get_front(numbers);
                auto array = enabledElementsArray<N>(tail);
                auto index = N - ctll::size(number);
                array[index] = number.value;
                return array;
            }
        }

        template<int ...Values>
        consteval auto toArray(ctll::list<Number<Values>...> numbers) noexcept
        {
            return toArrayHelper<ctll::size(numbers)>(numbers);
        }

        template<typename MaxFrequency, int AtomicNumber = 1>
        consteval auto enabledElements(auto smarts) noexcept
        {
            if constexpr (AtomicNumber > 104)
                return ctll::empty_list{};
            else if constexpr (smartsRequiresElement<MaxFrequency, AtomicNumber>(smarts))
                return ctll::push_front(Number<AtomicNumber>{}, enabledElements<MaxFrequency, AtomicNumber + 1>(smarts));
            else
                return enabledElements<MaxFrequency, AtomicNumber + 1>(smarts);
        }

        template<typename MaxFrequency>
        consteval auto enabledElementsArray(auto smarts) noexcept
        {
            return toArray(enabledElements<MaxFrequency>());
        }

        constexpr auto rejectElements(auto elements, Molecule::Molecule auto &mol) noexcept
        {
            if constexpr (ctll::empty(elements))
                return false;
            else {
                auto [element, tail] = ctll::pop_and_get_front(elements);
                for (auto atom : get_atoms(mol))
                    if (get_element(mol, atom) == element.value)
                        return rejectElements(tail, mol);

                return true;
            }
        }

        template<std::size_t N>
        constexpr auto rejectElements(const std::array<int, N> &elements, Molecule::Molecule auto &mol) noexcept
        {
            std::array<bool, N> found;
            std::ranges::fill(found, false);

            for (auto atom : get_atoms(mol))
                for (auto i = 0; i < N; ++i) {
                    if (found[i])
                        continue;
                    if (get_element(mol, atom) == elements[i])
                        found[i] = true;
                }

            return std::ranges::find(found, false) == found.end();
        }

    } // namespace impl

    template<typename MaxFrequency>
    struct ElementFilter
    {
        static consteval bool enable(auto smarts) noexcept
        {
            return !ctll::empty(impl::enabledElements<MaxFrequency>(smarts));
        }

        static constexpr bool reject(auto smarts, Molecule::Molecule auto &mol) noexcept
        {
            return impl::rejectElements(impl::enabledElements<MaxFrequency>(smarts), mol);
        }        
    };

} // namespace Kitimar::CTSmarts

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

    namespace impl {

        consteval auto enabledFilters(auto smarts, auto filters) noexcept
        {
            if constexpr (!ctll::empty(filters)) {
                auto [filter, tail] = ctll::pop_and_get_front(filters);
                if constexpr (filter.enable(smarts))
                    return ctll::push_front(filter, enabledFilters(smarts, tail));
                else
                    return enabledFilters(smarts, tail);
            } else
                return ctll::empty_list{};
        }

        template<typename ...Filter>
        constexpr bool rejectMolecule(auto smarts, Molecule::Molecule auto &mol, ctll::list<Filter...> filters) noexcept
        {
            return (Filter::reject(smarts, mol) || ...);
        }

    } // namespace impl

    template<typename SmartsT, typename ...Filter>
    struct FilterPolicyHelper
    {
        static constexpr inline auto filters = impl::enabledFilters(SmartsT{}, ctll::list<Filter...>{});

        static constexpr bool reject(Molecule::Molecule auto &mol) noexcept
        {
            return impl::rejectMolecule(SmartsT{}, mol, filters);
        }

        consteval FilterPolicyHelper() noexcept {}
        consteval FilterPolicyHelper(SmartsT, ctll::list<Filter...>) noexcept {}
    };

    template<typename ...Filter>
    struct FilterPolicy
    {
        static constexpr inline auto filters = ctll::list<Filter...>{};
    };

} // namespace Kitimar::CTSmarts

#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <variant>
#include <type_traits>
#include <unordered_set>
#include <cassert>

#define ISOMORPHISM_DEBUG 0

namespace Kitimar::CTSmarts {

    enum class SearchType
    {
        Single,
        Unique,
        All
    };

    template<SearchType T>
    using SearchTypeTag = std::integral_constant<SearchType, T>;

    static constexpr auto Single = SearchTypeTag<SearchType::Single>{};
    static constexpr auto Unique = SearchTypeTag<SearchType::Unique>{};
    static constexpr auto All    = SearchTypeTag<SearchType::All>{};

    template<ctll::fixed_string SMARTS>
    consteval bool requiresExplicitHydrogens() noexcept
    {
        using H = AliphaticAtom<1>;
        // Handle [!H]
        constexpr auto atoms = replaceExpr(Not<H>{}, AnyAtom{}, Smarts<SMARTS>::atoms);
        return containsExpr(H{}, atoms);
    }

    struct NoFilterPolicy : FilterPolicy<> {};

    struct UnconditionalFilterPolicy : FilterPolicy<NumAtomBondFilter> {};

    struct ConditionalFilterPolicy : FilterPolicy<> {};

    template<Molecule::Molecule Mol, typename SmartsT, SearchType SearchT, SeedType SeedT, typename Config = DefaultConfig>
    class Isomorphism
    {

        public:
            using Index = decltype(get_index(std::declval<Mol>(), get_atom(std::declval<Mol>(), 0)));
            using Map = IsomorphismMap<Index, SmartsT::numAtoms>;
            using Maps = IsomorphismMaps<Index, SmartsT::numAtoms>;

            static constexpr inline auto invalidIndex = static_cast<Index>(-1);
            static constexpr inline auto smarts = SmartsT{};
            static constexpr inline auto query = Config::Optimizer::create(smarts, SeedTypeTag<SeedT>{});

            static_assert(smarts.numBonds);
            static_assert(ctll::size(smarts.bonds) == ctll::size(query.bonds));

            constexpr const Map& map() const noexcept
            {
                return m_map.map;
            }

            // match

            bool match(const Mol &mol)
            {
                matchDfs(mol, nullptr);
                return isDone();
            }

            bool matchAtom(const Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                matchDfs(mol, nullptr, get_index(mol, atom));
                return isDone();
            }

            bool matchBond(const Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto sourceCallback = [this, targetIndex] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        setDone(true);
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                if (isDone() && m_map.contains(1, targetIndex))
                    return true;

                auto targetCallback = [this, sourceIndex] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        setDone(true);
                };
                matchDfs(mol, targetCallback, targetIndex);
                return isDone() && m_map.contains(1, sourceIndex);
            }

            // count

            auto count(const Mol &mol, int startAtom = -1)
            {
                auto n = 0;
                auto cb = [&n] (const auto &array) { ++n; };
                matchDfs(mol, cb, startAtom);
                return n;
            }

            auto countAtom(const Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                return count(mol, get_index(mol, atom));
            }

            auto countBond(const Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto n = 0;
                auto sourceCallback = [this, targetIndex, &n] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        ++n;
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                auto targetCallback = [this, sourceIndex, &n] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        ++n;
                };
                matchDfs(mol, targetCallback, targetIndex);
                return n;
            }

            // single

            auto single(const Mol &mol, int startAtom = -1)
            {
                matchDfs(mol, nullptr, startAtom);
                return std::make_tuple(isDone(), m_map.map());
            }

            auto singleAtom(const Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                return single(mol, get_index(mol, atom));
            }

            auto singleBond(const Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                auto sourceCallback = [this, targetIndex] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        setDone(true);
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                if (isDone() && m_map.contains(1, targetIndex))
                    return std::make_tuple(true, m_map.map());

                auto targetCallback = [this, sourceIndex] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        setDone(true);
                };
                matchDfs(mol, targetCallback, targetIndex);
                if (isDone() && m_map.contains(1, sourceIndex))
                    return std::make_tuple(true, m_map.map());

                return std::make_tuple(false, Map{});
            }

            // all

            Maps all(const Mol &mol, int startAtom = -1)
            {
                Maps maps;
                auto cb = [&maps] (const auto &map) {
                    maps.push_back(map);
                };
                matchDfs(mol, cb, startAtom);
                return maps;
            }

            auto allAtom(const Mol &mol, const auto &atom)
            {
                static_assert(SeedT == SeedType::Atom);
                return all(mol, get_index(mol, atom));
            }

            auto allBond(const Mol &mol, const auto &bond)
            {
                static_assert(SeedT == SeedType::Bond);
                reset(mol);
                auto source = get_source(mol, bond);
                auto target = get_target(mol, bond);
                auto sourceIndex = get_index(mol, source);
                auto targetIndex = get_index(mol, target);

                Maps maps;
                auto sourceCallback = [this, targetIndex, &maps] (const auto &map) {
                    if (m_map.contains(1, targetIndex))
                        maps.push_back(m_map.map());
                };
                matchDfs(mol, sourceCallback, sourceIndex);
                auto targetCallback = [this, sourceIndex, &maps] (const auto &map) {
                    if (m_map.contains(1, sourceIndex))
                        maps.push_back(m_map.map());
                };
                matchDfs(mol, targetCallback, targetIndex);
                return maps;
            }

        private:

            template<typename Arg, typename ...Args>
            void debug(Arg &&arg, Args &&...args)
            {
                #ifdef KITIMAR_WITH_IOSTREAM
                if constexpr (ISOMORPHISM_DEBUG) {
                    std::cout << std::forward<Arg>(arg);
                    ((std::cout << std::forward<Args>(args)), ...);
                    std::cout << '\n';
                }
                #endif // KITIMAR_WITH_IOSTREAM
            }

            bool matchAtom(const Mol &mol, const auto &atom, auto queryAtom) const noexcept
            {
                if (get_degree(mol, atom) < query.degrees[queryAtom.index])
                    return false;

                return matchAtomExpr(mol, atom, queryAtom.expr);
            }

            bool matchBond(const Mol &mol, const auto &bond, int queryBondIndex, auto queryBond) const noexcept
            {
                if constexpr (queryBond.isCyclic)
                    if (!is_ring_bond(mol, bond))
                        return false;

                return matchBondExpr(mol, bond, queryBond.expr);
            }

            void addAtom(int queryAtomIndex, auto atomIndex) noexcept
            {
                debug("    ", queryAtomIndex, " -> ", atomIndex);
                assert(static_cast<std::size_t>(queryAtomIndex) < m_map.map().size());
                assert(!m_map.containsQueryAtom(queryAtomIndex));
                m_map.add(queryAtomIndex, atomIndex);
            }

            void removeAtom(int queryAtomIndex, auto atomIndex) noexcept
            {
                debug("    backtrack: ", m_map(queryAtomIndex));
                assert(static_cast<std::size_t>(queryAtomIndex) < m_map.map().size());
                assert(m_map.containsQueryAtom(queryAtomIndex));
                m_map.remove(queryAtomIndex, atomIndex);
            }

            constexpr auto makeQueryBondInfo(auto queryBond) const noexcept
            {
                return std::make_tuple(queryBond.source, queryBond.target, queryBond.isCyclic, queryBond.isRingClosure);
            }

            constexpr auto getQueryBondInfo(int queryBondIndex) const noexcept
            {
                using R = std::tuple<int, int, bool, bool>;
                return with_n<ctll::size(query.bonds), R>(queryBondIndex, [this] (auto i) {
                    auto queryBond = get<i>(query.bonds);
                    return makeQueryBondInfo(queryBond);
                });
            }

            #ifdef KITIMAR_WITH_IOSTREAM
            void debugPoint(int queryBondIndex)
            {
                std::cout << "matchDfs<\"" << smarts.input() << "\">(queryBondIndex = " << queryBondIndex << "):" << std::endl;
                std::cout << "    m_map: " << m_map.map << std::endl;
                /*
                std::cout << "    m_mapped: [ ";
                for (auto v : m_mapped)
                    std::cout << v << " ";
                std::cout << "]" << std::endl;
                */

                if (queryBondIndex >= 0 && queryBondIndex < smarts.numBonds) {
                    // query
                    auto [querySource, queryTarget, isCyclic, isRingClosure] = getQueryBondInfo(queryBondIndex);
                    std::cout << "    query bond: " << querySource << " - " << queryTarget << std::endl;

                    auto source = m_map.map[querySource];
                    auto target = m_map.map[queryTarget];

                    std::cout << "    mol bond:   ";
                    if (source < 0)
                        std::cout << "? - ";
                    else
                        std::cout << source << " - ";
                    if (target < 0)
                        std::cout << "?" << std::endl;
                    else
                        std::cout << target << std::endl;
                }
                /*
                if (m_mapped.size() == m_map.size()) {
                    for (auto i = 0; i < m_mapped.size(); ++i)
                        assert(std::ranges::count(m_map, i) == m_mapped[i]);
                    assert(std::ranges::count(m_map, -1) == std::ranges::count(m_mapped, false));
                }
                */
            }
            #endif // KITIMAR_WITH_IOSTREAM

            template<typename Bonds = decltype(query.bonds)>
            void matchDfs(const Mol &mol, auto callback, int startAtom = -1, Bonds bonds = query.bonds)
            {
                if (isDone())
                    return;

                constexpr auto queryBondIndex = smarts.numBonds - ctll::size(bonds);
                if constexpr (ISOMORPHISM_DEBUG)
                    debugPoint(queryBondIndex);

                if constexpr (ctll::empty(bonds)) { // Found mapping?

                    debug("    found map: ", m_map.map());
                    addMapping(mol, callback);

                } else {
                    auto queryBond = ctll::front(bonds);
                    auto querySource = queryBond.source;
                    auto queryTarget = queryBond.target;

                    if constexpr (queryBond.isRingClosure) { // Ring closure?

                        auto source = get_atom(mol, m_map(querySource.index));
                        auto target = get_atom(mol, m_map(queryTarget.index));
                        if (!m_map.containsAtom(get_index(mol, target)))
                            return;
                        auto bond = Molecule::get_bond(mol, source, target);
                        if (bond == null_bond(mol))
                            return;
                        if (matchBondExpr(mol, bond, queryBond.expr))
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));

                    } else if (m_map.containsQueryAtom(querySource.index)) { // Source atom mapped?

                        auto source = get_atom(mol, m_map(querySource.index));

                        for (auto bond : get_bonds(mol, source)) {

                            auto target = Molecule::get_nbr(mol, bond, source);
                            auto targetIndex = get_index(mol, target);

                            if (m_map.containsAtom(targetIndex)) {
                                debug("    target already mapped: ", targetIndex);
                                continue;
                            }

                            // match bond
                            if (!matchBond(mol, bond, queryBondIndex, queryBond)) {
                                debug("    bond does not match");
                                continue;
                            }

                            // match target atom
                            if (!matchAtom(mol, target, queryTarget)) {
                                debug("    atom does not match");
                                continue;
                            }

                            // map target atom
                            addAtom(queryTarget.index, targetIndex);
                            matchDfs(mol, callback, startAtom, ctll::pop_front(bonds));
                            // exit as soon as possible if only one match is required
                            // (single mapping stored in m_map after returning)
                            if (isDone())
                                return;
                            // bracktrack target atom
                            if constexpr (!queryBond.isRingClosure)
                                removeAtom(queryTarget.index, targetIndex);

                        }

                    } else { // No mapped atoms

                        const auto numAtoms = num_atoms(mol);
                        if (numAtoms < smarts.numAtoms || num_bonds(mol) < smarts.numBonds)
                            return;

                        assert(!isDone());

                        assert(static_cast<std::size_t>(std::ranges::count(m_map.map(), invalidIndex)) == m_map.map().size());
                        //if constexpr (std::is_same_v<MappedPolicy<void>, MappedVector<void>>)
                        //    assert(std::ranges::count(m_mapped, true) == 0);

                        m_map.reset(numAtoms);

                        auto queryBond = ctll::front(query.bonds);
                        auto queryAtom = queryBond.source;

                        for (auto atom : get_atoms(mol)) {
                            if (startAtom != -1)
                                atom = get_atom(mol, startAtom);
                            auto index = get_index(mol, atom);

                            debug("    start atom: ",  index);

                            if (!matchAtom(mol, atom, queryAtom))
                                continue;

                            // map source atom, recursive dfs, backtrack
                            addAtom(queryAtom.index, index);
                            matchDfs(mol, callback, startAtom, query.bonds);
                            if (isDone())
                                return;
                            removeAtom(queryAtom.index, index);

                            if (startAtom != -1)
                                return;;
                        }

                    }
                }
            }

            constexpr auto reset(const Mol &mol) noexcept
            {
                //debug("reset()");
                setDone(false);
                m_map.reset(num_atoms(mol));
            }

            constexpr auto isDone() const noexcept
            {
                return m_done;
            }

            constexpr auto setDone(bool done) noexcept
            {
                m_done = done;
            }

            template<typename Callback>
            constexpr auto addMapping(const Mol &mol, Callback callback) noexcept
            {
                if constexpr (SearchT == SearchType::Single)
                    setDone(true);
                if constexpr (SearchT == SearchType::Unique) {
                    // create bit mask of atoms (to ensure uniqueness of mapping)
                    std::vector<bool> atoms(num_atoms(mol));
                    for (auto index : m_map.map())
                        atoms[index] = true;
                    // add the mapping to the result if it is unique
                    auto hash = std::hash<std::vector<bool>>()(atoms);
                    if (m_maps.find(hash) != m_maps.end())
                        return;
                    m_maps.insert(hash);
                }
                if constexpr (!std::is_same_v<std::nullptr_t, Callback>)
                    callback(m_map.map());
            }

            Config::template Map<Index, SmartsT::numAtoms> m_map; // current mapping: query atom index -> queried atom index
            std::conditional_t<SearchT == SearchType::Unique, std::unordered_set<std::size_t>, std::monostate> m_maps;
            bool m_done = false;
    };

    /*
    template<Molecule::Molecule Mol, typename SmartsT, MapType Type>
    class CentralAtomIsomorphism
    {
        public:
            using Map = IsomorphismMap<SmartsT::smarts>;
            using Maps = IsomorphismMaps<SmartsT::smarts>;

            static constexpr inline auto smarts = SmartsT{};

            static_assert(smarts.numBonds);

            CentralAtomIsomorphism()
            {
                m_map.fill(-1);
            }

            auto count(const Mol &mol, int startAtom = - 1)
            {
                auto n = 0;
                auto cb = [&n] (const auto &array) { ++n; };
                matchCentalAtom(mol, cb, startAtom);
                return n;
            }



        private:
            template<typename Bonds = decltype(smarts.bonds)>
            void matchCentralAtom(const Mol &mol, auto callback, int startAtom = -1, Bonds bonds = smarts.bonds)
            {
                auto queryBond = ctll::front(bonds);
                auto querySource = queryBond.source;
                auto queryTarget = queryBond.target;

                //if constexpr (queryTarget == smarts.centralAtom)

                for (auto atom : get_atoms(mol)) {
                    if (startAtom != -1)
                        atom = get_atom(mol, startAtom);
                    auto index = get_index(mol, atom);

                    debug("    start atom: ",  index);

                    if (matchAtom(mol, atom, querySource, queryBond.sourceExpr))
                        continue;

                    if (startAtom != -1)
                        return;;
                }
            }


            Map m_map; // current mapping: query atom index -> queried atom index
    };
    */

} // namespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::match<"SMARTS">(mol) -> bool

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr bool match(const Mol &mol)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol)) // FIXME: use std::ranges::find_if -> check assembly?
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return true;
            return false;
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol))
                if (impl::singleBondMatch(smarts, mol, bond))
                    return true;
            return false;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Molecule, Config>{};
            return iso.match(mol);
        }
    }

    //
    // Atom
    //

    // ctse::match_atom<"SMARTS">(mol, atom) -> bool

    namespace impl {

        template<typename SmartsT, typename Config, Molecule::Molecule Mol>
        constexpr bool match_atom(const Mol &mol, const auto &atom)
        {
            auto smarts = Config::transformSmarts(SmartsT{}); // FIXME: replaceExpr(BasicSmarts<...>)
            if constexpr (smarts.isSingleAtom) {
                // Optimize single atom SMARTS
                return impl::singleAtomMatch(smarts, mol, atom);
            } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
                // Optimize single bond SMARTS
                if (!matchAtomExpr(mol, atom, get<0>(smarts.atoms).expr))
                    return false;
                for (auto bond : get_bonds(mol, atom)) {
                    if (!matchBondExpr(mol, bond, get<0>(smarts.bonds).expr))
                        continue;
                    if (matchAtomExpr(mol, Kitimar::Molecule::get_nbr(mol, bond, atom), get<1>(smarts.atoms).expr))
                        return true;
                }
                return false;
            } else {
                auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Atom, Config>{};
                return iso.matchAtom(mol, atom);
            }
        }

    } // namespace impl

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr bool match_atom(const Mol &mol, const auto &atom)
    {
        return impl::match_atom<Smarts<SMARTS>, Config, Mol>(mol, atom);

    }

    //
    // Bond
    //

    // ctse::match_bond<"SMARTS">(mol, bond) -> bool

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr bool match_bond(const Mol &mol, const auto &bond)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            return impl::singleBondMatch(smarts, mol, bond);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::All, SeedType::Bond, Config>{};
            return iso.matchBond(mol, bond);
        }
    }

    //
    // Atom/Bond
    //

    // ctse::match<"SMARTS">(mol, atom) -> bool

    CTSMARTS_API_OVERLOAD_ATOM(match)

    // ctse::match<"SMARTS">(mol, bond) -> bool

    CTSMARTS_API_OVERLOAD_BOND(match)

} // mamespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::count<"SMARTS">(mol, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    constexpr auto count(const Mol &mol, SearchTypeTag<M> = {})
    {
        static_assert(M != SearchType::Single, "Use CTSmarts::contains<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            auto n = 0;
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    ++n;
            return n;
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            auto n = 0;
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == SearchType::Unique) {
                    if (impl::singleBondMatch(smarts, mol, bond))
                        ++n;
                } else {
                    n += impl::singleBondCount(smarts, mol, bond);
                }
            }
            return n;
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Molecule, Config>{};
            return iso.count(mol);
        }
    }

    // ctse::count_unique<"SMARTS">(mol) -> std::integeral

    CTSMARTS_API_UNIQUE(count)

    // ctse::count_all<"SMARTS">(mol) -> std::integeral

    CTSMARTS_API_ALL(count)

    //
    // Atom
    //

    // ctse::count_atom<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    constexpr auto count_atom(const Mol &mol, const auto &atom, SearchTypeTag<M> = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        if constexpr (smarts.isSingleAtom)
            return match_atom<SMARTS, Config>(mol, atom) ? 1 : 0;
        else {
            auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Atom, Config>{};
            return iso.countAtom(mol, atom);
        }
    }

    // ctse::count_atom_unique<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_ATOM_UNIQUE(count)

    // ctse::count_atom_all<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_ATOM_ALL(count)

    //
    // Bond
    //

    // ctse::count_bond<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::integeral

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    constexpr auto count_bond(const Mol &mol, const auto &bond, SearchTypeTag<M> = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            if constexpr (M == SearchType::Unique) {
                if (impl::singleBondMatch(smarts, mol, bond))
                    return 1;
            } else {
                return impl::singleBondCount(smarts, mol, bond);
            }
            return 0;
        }
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Bond, Config>{};
        return iso.countBond(mol, bond);
    }

    // ctse::count_bond_unique<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_BOND_UNIQUE(count)

    // ctse::count_bond_all<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_BOND_ALL(count)

    //
    // Atom/Bond
    //

    // ctse::count<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::integeral

    CTSMARTS_API_OVERLOAD_ATOM_SEARCH(count)

    // ctse::count<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::integeral

    CTSMARTS_API_OVERLOAD_BOND_SEARCH(count)

    // ctse::count_unique<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(count)

    // ctse::count_unique<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_OVERLOAD_BOND_UNIQUE(count)

    // ctse::count_all<"SMARTS">(mol, atom) -> std::integeral

    CTSMARTS_API_OVERLOAD_ATOM_ALL(count)

    // ctse::count_all<"SMARTS">(mol, bond) -> std::integeral

    CTSMARTS_API_OVERLOAD_BOND_ALL(count)

} // mamespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::map<"SMARTS">(mol) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto map(const Mol &mol)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        using Map = IsomorphismMap<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, Map{get_index(mol, atom)});
            return std::make_tuple(false, Map{});
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                auto m = impl::singleBondMap<Map>(smarts, mol, bond);
                if (std::get<0>(m))
                    return m;
            }
            return std::make_tuple(false, Map{});
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Molecule, Config>{};
            return iso.single(mol);
        }
    }

    //
    // Atom
    //

    // ctse::map_atom<"SMARTS">(mol, atom) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto map_atom(const Mol &mol, const auto &atom)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Atom, Config>{};
        return iso.singleAtom(mol, atom);
    }

    //
    // Bond
    //

    // ctse::map_bond<"SMARTS">(mol, bond) -> std::tuple<bool, std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    constexpr auto map_bond(const Mol &mol, const auto &bond)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(smarts.numBonds, "There should at least be one bond in the SMARTS expression.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Bond, Config>{};
        return iso.singleBond(mol, bond);
    }

    //
    // Atom/Bond
    //

    // ctse::map<"SMARTS">(mol, atom) -> std::tuple<bool, std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM(map)

    // ctse::map<"SMARTS">(mol, bond) -> std::tuple<bool, std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND(map)

} // mamespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::maps<"SMARTS">(mol, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    constexpr auto maps(const Mol &mol, SearchTypeTag<M> = {})
    {
        static_assert(M != SearchType::Single, "Use CTSmarts::map<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        using Maps = IsomorphismMaps<decltype(get_index(mol, get_atom(mol, 0))), smarts.numAtoms>;
        using Map = Maps::value_type;

        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            Maps maps;
            maps.reserve(num_atoms(mol));
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    maps.push_back(Map{get_index(mol, atom)});
            return maps;
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            Maps maps;
            maps.reserve(num_bonds(mol));
            for (auto bond : get_bonds(mol)) {
                if constexpr (M == SearchType::Unique) {
                    auto [found, map] = impl::singleBondMap<Map>(smarts, mol, bond);
                    if (found)
                        maps.push_back(map);
                } else {
                    impl::singleBondMaps(smarts, mol, bond, maps);
                }
            }
            return maps;
        //} else if constexpr (smarts.centralAtom != -1) {

        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Molecule, Config>{};
            return iso.all(mol);
        }
    }

    // ctse::maps_unique<"SMARTS">(mol) -> std::vector<std::array<int, N>>

    CTSMARTS_API_UNIQUE(maps)

    // ctse::maps_all<"SMARTS">(mol) -> std::vector<std::array<int, N>>

    CTSMARTS_API_ALL(maps)

    //
    // Atom
    //

    // ctse::maps_atom<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    constexpr auto maps_atom(const Mol &mol, const auto &atom, SearchTypeTag<M> SearchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::map_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Atom, Config>{};
        return iso.allAtom(mol, atom);
    }

    // ctse::maps_atom_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_ATOM_UNIQUE(maps)

    // ctse::maps_atom_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_ATOM_ALL(maps)

    //
    // Bond
    //

    // ctse::maps_bond<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    template<ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    constexpr auto maps_bond(const Mol &mol, const auto &bond, SearchTypeTag<M> SearchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom /*&& !smarts.isSingleBond*/,
                "Use CTSmarts::map_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Bond, Config>{};
        return iso.allBond(mol, bond);
    }

    // ctse::maps_bond_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_BOND_UNIQUE(maps)

    // ctse::maps_bond_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_BOND_ALL(maps)

    //
    // Atom/Bond
    //

    // ctse::maps<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_SEARCH(maps)

    // ctse::maps<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_SEARCH(maps)

    // ctse::maps_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(maps)

    // ctse::maps_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_UNIQUE(maps)

    // ctse::maps_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_ALL(maps)

    // ctse::maps_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_ALL(maps)

} // mamespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::capture<"SMARTS">(mol) -> std::tuple<bool, Atom...>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto capture(const Mol &mol)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    return std::make_tuple(true, atom);
            return std::make_tuple(false, null_atom(mol));
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            for (auto bond : get_bonds(mol)) {
                auto matchType = impl::singleBondMatch(smarts, mol, bond);
                if (matchType)
                    return impl::singleBondCapture(smarts, mol, bond, matchType);
            }
            return impl::singleBondCapture(smarts, mol, null_bond(mol), 0);
        } else {
            auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Molecule, Config>{};
            constexpr auto captureSet = captureMapping(smarts);
            auto [found, map] = iso.single(mol);
            return impl::captureMatchAtoms(mol, smarts, captureSet, found, map);
        }
    }

    //
    // Atom
    //

    // ctse::capture_atom<"SMARTS">(mol, atom) -> std::tuple<bool, Atom...>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto capture_atom(const Mol &mol, const auto &atom)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Atom, Config>{};
        constexpr auto captureSet = captureMapping(smarts);
        auto [found, map] = iso.singleAtom(mol, atom);
        return impl::captureMatchAtoms(mol, smarts, captureSet, found, map);
    }

    //
    // Bond
    //

    // ctse::capture_bond<"SMARTS">(mol, bond) -> std::tuple<bool, Atom...>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, Molecule::Molecule Mol>
    auto capture_bond(const Mol &mol, const auto &atom)
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(!smarts.isSingleAtom, "Use CTSmarts::match_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), SearchType::Single, SeedType::Bond, Config>{};
        constexpr auto captureSet = captureMapping(smarts);
        auto [found, map] = iso.singleBond(mol, atom);
        return impl::captureMatchAtoms(mol, smarts, captureSet, found, map);
    }

    //
    // Atom/Bond
    //

    // ctse::capture<"SMARTS">(mol, atom) -> std::tuple<bool, Atom...>

    CTSMARTS_API_OVERLOAD_ATOM(capture)

    // ctse::capture<"SMARTS">(mol, bond) -> std::tuple<bool, Atom...>

    CTSMARTS_API_OVERLOAD_BOND(capture)

} // mamespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::captures<"SMARTS">(mol, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    auto captures(const Mol &mol, SearchTypeTag<M> = {})
    {
        static_assert(M != SearchType::Single, "Use CTSmarts::capture<\"SMARTS\">(mol) to check for a single match.");
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static constexpr auto captureSet = captureMapping(smarts);
        using IndexMap = std::array<decltype(get_index(mol, get_atom(mol, 0))), captureSet.size() ? captureSet.size() : smarts.numAtoms>;
        using AtomMap = std::array<decltype(get_atom(mol, 0)), captureSet.size() ? captureSet.size() : smarts.numAtoms>;
        using AtomMaps = std::vector<AtomMap>;

        if constexpr (smarts.isSingleAtom) {
            // Optimize single atom SMARTS
            AtomMaps maps;
            maps.reserve(num_atoms(mol));
            for (auto atom : get_atoms(mol))
                if (impl::singleAtomMatch(smarts, mol, atom))
                    maps.push_back(impl::toCapture(mol, smarts, captureSet, true, IndexMap{get_index(mol, atom)}));
            return maps;
        } else if constexpr ((Config::specialize & Specialize::Bond) && smarts.isSingleBond) {
            // Optimize single bond SMARTS
            AtomMaps maps;
            maps.reserve(num_bonds(mol));
            if constexpr (captureSet.size() == smarts.numAtoms) {
                for (auto bond : get_bonds(mol))
                    impl::singleBondCaptures(smarts, mol, bond, captureSet, maps, M == SearchType::Unique);
            } else {
                for (auto bond : get_bonds(mol))
                    impl::singleBondCaptures(smarts, mol, bond, captureSet, maps, false);

                // Remove duplicates
                if constexpr (M == SearchType::Unique) {
                    auto proj = [&mol] (const auto &atoms) { return impl::captureHash(mol, atoms); };
                    std::ranges::sort(maps, {}, proj);
                    const auto [first, last] = std::ranges::unique(maps, {}, proj);
                    maps.erase(first, last);
                }
            }
            return maps;
        } else {
            if constexpr (captureSet.size() == smarts.numAtoms) {
                auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Molecule, Config>{};
                if constexpr (__cpp_lib_ranges >= 202110L)
                    return iso.all(mol) | std::views::transform([&] (const auto &map) {
                        return impl::toCapture(mol, smarts, true, map, captureSet);
                    });
                else
                    // missing std::ranges::owning_view
                    return impl::toCaptures(mol, iso, captureSet, iso.all(mol));
            } else {
                auto iso = Isomorphism<Mol, decltype(smarts), SearchType::All, SeedType::Molecule, Config>{};
                AtomMaps maps = impl::toCaptures(mol, iso, captureSet, iso.all(mol));

                // Remove duplicates
                if constexpr (M == SearchType::Unique) {
                    auto proj = [&mol] (const auto &atoms) { return impl::captureHash(mol, atoms); };
                    std::ranges::sort(maps, {}, proj);
                    const auto [first, last] = std::ranges::unique(maps, {}, proj);
                    maps.erase(first, last);
                }

                return maps;
            }
        }
    }

    // ctse::captures_unique<"SMARTS">(mol) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_UNIQUE(captures)

    // ctse::captures_all<"SMARTS">(mol) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_ALL(captures)

    //
    // Atom
    //

    // ctse::captures_atom<"SMARTS">(mol, atom, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    auto captures_atom(const Mol &mol, const auto &atom, SearchTypeTag<M> searchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::capture_atom<\"SMARTS\">(mol, atom) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Atom, Config>{};
        static constexpr auto captureSet = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.allAtom(mol, atom) | std::views::transform([&] (const auto &map) {
                return impl::toCapture(mol, smarts, captureSet, true, map);
            });
        else
            // missing std::ranges::owning_view
            return impl::toCaptures(mol, iso, captureSet, iso.allAtom(mol, atom));
    }

    // ctse::captures_atom_unique<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_ATOM_UNIQUE(captures)

    // ctse::captures_atom_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_ATOM_ALL(captures)

    //
    // Bond
    //

    // ctse::captures_bond<"SMARTS">(mol, bond, CTSmarts::[Unique, All]) -> std::vector<std::array<Atom, N>>

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig, SearchType M = SearchType::Unique, Molecule::Molecule Mol>
    auto captures_bond(const Mol &mol, const auto &bond, SearchTypeTag<M> searchType = {})
    {
        auto smarts = Config::transformSmarts(Smarts<SMARTS>{});
        static_assert(M != SearchType::Single && !smarts.isSingleAtom,
                    "Use CTSmarts::capture_bond<\"SMARTS\">(mol, bond) to check for a single match.");
        auto iso = Isomorphism<Mol, decltype(smarts), M, SeedType::Bond, Config>{};
        static constexpr auto captureSet = captureMapping(smarts);
        if constexpr (__cpp_lib_ranges >= 202110L)
            return iso.allBond(mol, bond) | std::views::transform([&] (const auto &map) {
                return impl::toCapture(mol, smarts, captureSet, true, map);
            });
        else
            // missing std::ranges::owning_view
            return impl::toCaptures(mol, iso, captureSet, iso.allBond(mol, bond));
    }

    // ctse::captures_bond_unique<"SMARTS">(mol, bond) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_BOND_UNIQUE(captures)

    // ctse::captures_bond_all<"SMARTS">(mol, atom) -> std::vector<std::array<Atom, N>>

    CTSMARTS_API_BOND_ALL(captures)

    //
    // Atom/Bond
    //

    // ctse::captures<"SMARTS">(mol, atom, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_SEARCH(captures)

    // ctse::captures<"SMARTS">(mol, bond, ctse::[Unique, All]) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_SEARCH(captures)

    // ctse::captures_unique<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_UNIQUE(captures)

    // ctse::captures_unique<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_UNIQUE(captures)

    // ctse::captures_all<"SMARTS">(mol, atom) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_ATOM_ALL(captures)

    // ctse::captures_all<"SMARTS">(mol, bond) -> std::vector<std::array<int, N>>

    CTSMARTS_API_OVERLOAD_BOND_ALL(captures)

} // mamespace Kitimar::CTSmarts

namespace Kitimar::CTSmarts {

    //
    // Molecule
    //

    // ctse::find<"SMARTS">(mol) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find(const Molecule::Molecule auto &mol)
    {
        auto caps = capture<SMARTS, Config>(mol);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    // ctse::find_unique<"SMARTS">(mol) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find_unique(const Molecule::Molecule auto &mol)
    {
        return captures_unique<SMARTS, Config>(mol) | std::views::transform([] (const auto &capture) {
            return capture[0];
        });
    }

    //
    // Atom
    //

    // ctse::find_atom<"SMARTS">(mol, atom) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find_atom(const Molecule::Molecule auto &mol, const auto &atom)
    {
        auto caps = capture_atom<SMARTS, Config>(mol, atom);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // Bond
    //

    // ctse::find_bond<"SMARTS">(mol, bond) -> Atom (null atom if there is no match)

    template <ctll::fixed_string SMARTS, typename Config = DefaultConfig>
    auto find_bond(const Molecule::Molecule auto &mol, const auto &bond)
    {
        auto caps = capture_bond<SMARTS, Config>(mol, bond);
        static_assert(std::tuple_size<decltype(caps)>{} >= 2); // FIXME: better validation 2 or # SMARTS atoms
        return std::get<1>(caps);
    }

    //
    // Atom/Bond
    //

    // ctse::find_atom<"SMARTS">(mol, atom) -> Atom (null atom if there is no match)

    CTSMARTS_API_OVERLOAD_ATOM(find)

    // ctse::find_bond<"SMARTS">(mol, bond) -> Atom (null atom if there is no match)

    CTSMARTS_API_OVERLOAD_BOND(find)

} // mamespace Kitimar::CTSmarts

/*

API
===

Match
-----

    match(mol)
    match_atom(mol, atom)
    match_bond(mol, bond)
    match(mol, atom/bond)

    Return type: bool

Count
-----

    count(mol, type = Unique)
    count_unique(mol)
    count_all(mol)

    count_atom(mol, atom, type = Unique)
    count_atom_unique(mol, atom)
    count_atom_all(mol, atom)

    count_bond(mol, bond, type = Unique)
    count_bond_unique_bond(mol, bond)
    count_bond_all(mol, bond)

    count(mol, atom/bond, type = Unique)
    count_unique(mol, atom/bond)
    count_all(mol, atom/bond)

    Return type: int

Map
---

    map(mol)
    map_atom(mol, atom)
    map_bond(mol, bond)
    map(mol, atom/bond)

    Return type: std::tuple<bool, std::array<int, N>>

    Example
    ~~~~~~~

        auto [found, map] = ctse::map<"C=O">(mol);
        if (found) {
            // Use map...
        }

Maps
----

    maps(mol, type = Unique)
    maps_unique(mol)
    maps_all(mol)

    maps_atom(mol, atom, type = Unqiue)
    maps_atom_unique(mol, atom)
    maps_atom_all(mol, atom)

    maps_bond(mol, bond, type = Unique)
    maps_bond_unique(mol, bond)
    maps_bond_all(mol, bond)

    maps(mol, atom/bond, type = Unique)
    maps_unique(mol, atom/bond)
    maps_all(mol, atom/bond)

    Return type: std::vector<std::array<int, N>>

Capture
-------

    capture(mol)
    capture_atom(mol, atom)
    capture_bond(mol, bond)
    capture(mol, atom/bond)

    Return type: std::tuple<bool, Atom...>

    Examples
    ~~~~~~~~

        auto [found, C, O] = ctse::capture<"C=O">
        if (found) {
            // Use atoms C and O...
        }

        auto [found, C, O] = ctse::capture<"N[C:1]=[O:2]">
        if (found) {
            // Use atoms C and O...
        }

Captures
--------

    captures(mol, type = Unique)
    captures_unique(mol)
    captures_all(mol)

    captures_atom(mol, atom, type = Unique)
    captures_atom_unique(mol, atom)
    captures_atom_all(mol, atom)

    captures_bond(mol, bond, type = Unique)
    captures_bond_unique(mol)
    captures_bond_all(mol)

    captures(mol, atom/bond, type = Unique)
    captures_unique(mol, atom/bond)
    captures_all(mol, atom/bond)

    Return type: std::vector<std::tuple<Atom...>>

    Example
    ~~~~~~~

        for (auto [C, N] : ctse::captures<"C-N">(mol)) {
            // Use atoms C and N
        }

*/

namespace ctse = Kitimar::CTSmarts;
