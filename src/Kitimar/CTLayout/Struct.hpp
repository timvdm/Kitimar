#pragma once

#include <Kitimar/CTLayout/Source.hpp>

#include <ctll/list.hpp>

#include <type_traits>
#include <cstdint>
#include <concepts>
#include <cassert>
#include <array>
#include <algorithm>

namespace Kitimar::CTLayout {

    template<typename T, typename U, typename ...Us>
    constexpr auto indexOf(T, ctll::list<U, Us...>) noexcept
    {
        if constexpr (std::is_same_v<T, U>)
            return 0;
        else
            return 1 + indexOf(T{}, ctll::list<Us...>{});
    }


    // isFixedSize(ctll::list)

    constexpr bool isFixedSize(ctll::empty_list) noexcept { return true; }

    template<typename T, typename ...Ts>
    constexpr bool isFixedSize(ctll::list<T, Ts...>) noexcept
    {
        if constexpr (!isFixedSize(T{}))
            return false;
        else
            return isFixedSize(ctll::list<Ts...>{});
    }

    // contains(ctll::list)

    template<typename U>
    constexpr bool contains(ctll::empty_list, U) noexcept { return false; }

    template<typename U, typename T, typename ...Ts>
    constexpr bool contains(ctll::list<T, Ts...>, U) noexcept
    {
        if constexpr (contains(T{}, U{}))
            return true;
        return contains(ctll::list<Ts...>{}, U{});
    }

    // sizeOf(ctll::list)

    template<typename T, typename ...Ts, typename Source = PtrSource>
    constexpr std::size_t sizeOf(ctll::list<T, Ts...>, Source data = {}) noexcept
    {
        auto sizeOfT = sizeOf(T{}, data);
        if constexpr(sizeof...(Ts))
            if constexpr (isFixedSize(ctll::list<Ts...>{}))
                return sizeOfT + sizeOf(ctll::list<Ts...>{});
            else
                return sizeOfT + sizeOf(ctll::list<Ts...>{}, data + sizeOfT);
        else
            return sizeOfT;
    }

    // offset(ctll::list) - FixedSize

    constexpr std::size_t offset(ctll::empty_list, auto) noexcept { return 0; }

    template<typename U, typename T, typename ...Ts>
    constexpr std::size_t offset(ctll::list<T, Ts...>, U) noexcept
    {
        if constexpr (contains(T{}, U{}))
            return offset(T{}, U{});
        else
            return sizeOf(T{}) + offset(ctll::list<Ts...>{}, U{});
    }

    // offset(ctll::list) - VariableSize

    constexpr std::size_t offset(ctll::empty_list, auto, auto,  auto) noexcept { return 0; }

    template<typename U, typename T, typename ...Ts>
    constexpr std::size_t offset(ctll::list<T, Ts...>, U, auto data, auto index) noexcept
    {
        if constexpr (contains(T{}, U{})) {
            if constexpr (isFixedSize(T{}))
                return offset(T{}, U{});
            else
                return offset(T{}, U{}, data, index);
        } else if constexpr (isFixedSize(T{})) {
            auto sizeOfT = sizeOf(T{});
            return sizeOfT + offset(ctll::list<Ts...>{}, U{}, data + sizeOfT, index);
        } else {
            auto sizeOfT = sizeOf(T{}, data);
            return sizeOfT + offset(ctll::list<Ts...>{}, U{}, data + sizeOfT, index);
        }


    }

    // FXIME: Ts = Value, Struct, Array
    template<typename ...Ts>
    struct Struct
    {
        using Members = ctll::list<Ts...>;
        static constexpr inline auto members = Members{}; // FIXME: remove?

        static constexpr std::size_t size() noexcept { return sizeOf(Members{}, nullptr); } // FIXME: remove

        constexpr Struct() noexcept {}
        constexpr Struct(Members) noexcept {}
    };

    template<typename T>
    concept IsStruct = std::derived_from<T, decltype(Struct{T::members})>;

    // isStruct

    template<typename ...Ts>
    constexpr auto isStruct(Struct<Ts...>) { return true; }
    constexpr auto isStruct(...) { return false; }

    // isFixedSize

    template<typename ...Ts>
    constexpr bool isFixedSize(Struct<Ts...>) noexcept
    {
        return isFixedSize(ctll::list<Ts...>{});
    }

    // contains

    template<typename StructT, typename T>
    constexpr bool contains(StructT, T) noexcept requires IsStruct<StructT>
    {
        if constexpr (std::is_same_v<StructT, T>)
            return true;
        return contains(StructT::members, T{});
    }

    // sizeOf

    template<typename ...Ts, typename Source = PtrSource>
    constexpr std::size_t sizeOf(Struct<Ts...>, Source data = {}) noexcept
    {
        if constexpr (!isFixedSize(Struct<Ts...>{}))
            assert(data);
        return /*sizeof(std::size_t) +*/ sizeOf(ctll::list<Ts...>{}, data); // FIXME: optimize? -> yes, quick skip
    }

    // offset

    template<typename StructT, typename T, typename Source = PtrSource, std::integral Index = uint32_t>
    constexpr std::size_t offset(StructT, T, Source data = {}, Index index = {}) noexcept requires IsStruct<StructT>
    {
        if constexpr (std::is_same_v<StructT, T>)
            return 0;
        if constexpr (isFixedSize(StructT{})) {
            return offset(StructT::members, T{});
        } else {
            assert(data);
            if constexpr (!contains(StructT{}, T{}))
                return sizeOf(StructT{}, data);
            else
                return offset(StructT::members, T{}, data, index);
        }
    }

    // StructObject

    template<typename StructT, typename SourceT = PtrSource>
    class StructObject
    {
        public:
            using Type = StructT;
            using Members = typename StructT::Members;
            using Source = SourceT;

            constexpr StructObject(StructT, Source data) noexcept : m_data(data)
            {
            }

            template<typename T>
            constexpr auto get(T) const noexcept
            {
                static_assert(contains(StructT{}, T{}));
                auto off = offset(StructT{}, T{}, m_data);
                return toObject(T{}, m_data + off);
            }

            constexpr auto data() const noexcept
            {
                return m_data;
            }

        private:
            Source m_data = {};
    };

    // toObject

    template<typename StructT, typename Source = PtrSource>
    constexpr auto toObject(StructT, Source data) noexcept requires IsStruct<StructT>
    {
        return StructObject{StructT{}, data};
    }

    // StructWriter

    template<typename StructT, typename Sink, typename Finished = std::nullptr_t>
    class StructWriter
    {
        public:
            using Type = StructT;
            using Members = typename StructT::Members;

            StructWriter(StructT, Sink sink, Finished finished = nullptr) : m_sink(sink), m_finished{finished}
            {
            }

            template<typename T>
            constexpr auto get(T) noexcept
            {
                static_assert(contains(StructT{}, T{}));
                constexpr auto index = indexOf(T{}, Members{});
                assert(!m_locked[index]);
                m_locked[index] = true;
                auto off = offset(Type{}, T{}, m_sink);
                return toWriter(T{}, m_sink + off, [this] () {
                    m_written[indexOf(T{}, Members{})] = true;
                    if constexpr (!std::is_same_v<Finished, std::nullptr_t>)
                        if (std::ranges::count(m_written, true) == m_written.size())
                            m_finished();
                });
            }

        private:
            Sink m_sink;
            std::array<bool, ctll::size(Type::members)> m_written = {};
            std::array<bool, ctll::size(Type::members)> m_locked = {};
            Finished m_finished;

    };

    // toWriter

    template<typename StructT, typename Sink, typename Finished = std::nullptr_t>
    constexpr auto toWriter(StructT, Sink sink, Finished finished = nullptr) noexcept requires IsStruct<StructT>
    {
        return StructWriter{StructT{}, sink, finished};
    }

} // namespace Kitimar::CTLayout
