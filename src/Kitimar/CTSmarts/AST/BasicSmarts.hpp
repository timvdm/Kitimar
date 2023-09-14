#pragma once

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
