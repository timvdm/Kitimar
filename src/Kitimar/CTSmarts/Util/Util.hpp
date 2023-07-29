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
    struct Number : std::integral_constant<int, Value> {};

    template<int Value>
    struct Char : std::integral_constant<char, Value> {};

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
