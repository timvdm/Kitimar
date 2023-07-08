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
    struct Number
    {
        static constexpr inline auto value = Value;
    };

    template<int Value>
    struct Char
    {
        static constexpr inline auto value = Value;
    };

    //
    // ctll::list extensions
    //

    template<int Index, typename ...Ts>
    constexpr auto get(ctll::list<Ts...>)
    {
        return std::get<Index>(std::tuple<Ts...>());
    }

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


    template<int ...N>
    constexpr auto toArray(ctll::list<Number<N>...>)
    {
        return std::array<int, sizeof...(N)>({ N... });
    }

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
