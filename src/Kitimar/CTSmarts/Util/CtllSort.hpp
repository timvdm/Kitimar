#pragma once

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
            if constexpr (Compare{}(Project{}(decltype(a){}), Project{}(decltype(b){})))
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
