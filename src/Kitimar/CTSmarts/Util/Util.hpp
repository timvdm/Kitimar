#pragma once

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

    namespace detail {

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
            return std::array{ detail::with_integral_constant<I, R, F&&>... };
        }(std::make_index_sequence<N>{});

        return invokeArray.at(n)(std::forward<F>(f));
    }



} // namespace ctsmarts
