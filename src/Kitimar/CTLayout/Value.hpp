#pragma once

#include "Source.hpp"

#include <type_traits>
#include <concepts>

namespace Kitimar::CTLayout {

    // FXIME: T = integral or floating_point
    template<typename T>
    struct Value
    {
        using Type = T;

        static constexpr std::size_t size() noexcept { return sizeof(T); } // FIXME: remove...

    };

    template<typename T>
    concept IsValue = std::derived_from<T, Value<typename T::Type>>;

    // isValue

    template<typename T>
    constexpr bool isValue(Value<T>) noexcept { return true; }
    constexpr bool isValue(...) noexcept { return false; }

    // isFixedSize

    template<typename T>
    constexpr bool isFixedSize(Value<T>) noexcept { return true; }

    // contains

    template<typename ValueT, typename T>
    constexpr bool contains(ValueT, T) noexcept requires IsValue<ValueT>
    {
        return std::is_same_v<ValueT, T>;
    }

    // sizeOf

    template<typename T, typename Source = PtrSource>
    constexpr std::size_t sizeOf(Value<T>, Source = {}) noexcept { return sizeof(T); }

    // offset

    template<typename ValueT, typename T>
    constexpr std::size_t offset(ValueT, T) noexcept requires IsValue<ValueT>
    {
        if constexpr (contains(ValueT{}, T{}))
            return 0;
        return sizeOf(ValueT{});
    }

    // ValueObject

    template<typename ValueT, typename Source = PtrSource>
    class ValueObject
    {
        public:
            using Type = typename ValueT::Type;

            constexpr ValueObject(ValueT, Source data) noexcept : m_data(data)
            {
            }

            constexpr Type get() const noexcept
            {
                return m_data.template read<Type>();
            }

            constexpr void set(Type value) noexcept
            {
                m_data.write(value);
            }

        private:
            Source m_data = {};
    };

    // toObject

    template<typename ValueT, typename Source = PtrSource>
    constexpr auto toObject(ValueT, Source data) noexcept requires IsValue<ValueT>
    {
        return ValueObject{ValueT{}, data};
    }

    // ValueWriter

    template<typename ValueT, typename Sink, typename Finished = std::nullptr_t>
    class ValueWriter
    {
        public:
            using Type = typename ValueT::Type;

            ValueWriter(ValueT, Sink sink, Finished finished = nullptr) : m_sink(sink), m_finished{finished}
            {
            }

            constexpr void write(Type value) noexcept
            {
                m_sink.write(value);
                m_written = true;
                if constexpr (!std::is_same_v<Finished, std::nullptr_t>)
                    m_finished();
            }

        private:
            Sink m_sink;
            Finished m_finished;
            bool m_written;
    };

    // toWriter

    template<typename ValueT, typename Sink, typename Finished = std::nullptr_t>
    constexpr auto toWriter(ValueT, Sink sink, Finished finished = nullptr) noexcept requires IsValue<ValueT>
    {
        return ValueWriter{ValueT{}, sink, finished};
    }

} // namespace Kitimar::CTLayout
