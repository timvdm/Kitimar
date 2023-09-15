#pragma once

#include "../Optimizer/ExpressionFrequency.hpp"
#include "../Util/Ctll.hpp"

#include "../../Molecule/Molecule.hpp"

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
