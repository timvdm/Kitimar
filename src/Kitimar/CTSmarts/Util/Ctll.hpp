#pragma once

#include "Util.hpp"

#include <ctll/list.hpp>

namespace Kitimar::CTSmarts {

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

} // namespace Kitimar::CTSmarts
