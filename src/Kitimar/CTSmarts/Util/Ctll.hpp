#pragma once

#include "Util.hpp"

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
