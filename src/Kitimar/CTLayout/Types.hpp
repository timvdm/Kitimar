#pragma once

#include <Kitimar/CTLayout/Value.hpp>
#include <Kitimar/CTLayout/Struct.hpp>
#include <Kitimar/CTLayout/Vector.hpp>

#include <ctll/list.hpp>

#include <cstdint>
#include <utility>

/*
 * - Only Layout can contain arrays
 *   - Value and Struct have fixed size
 *
 *
 * - Value
 * - Struct
 *   - Value
 *   - Struct
 * - Array
 *   - Value
 *   - Struct
 * - Layout
 *   - Value
 *   - Struct
 *   - Array
 *
 *
 *
 * sizeOf
 * ------
 *
 * Value<T>         sizeof(T)
 * Struct<Ts...>    sizeOf(Ts...)
 * Array<T>
 *
 *
 *
 *
 *
 */


namespace Kitimar::CTLayout {

































    using SizeT = uint32_t;




    /*
    // FXIME: T = integral or floating_point
    template<typename T>
    struct Value
    {
        using Type = T;
        static constexpr inline auto type = Type{};

        static constexpr std::size_t size() noexcept { return sizeof(Type); }
    };
    */

    /*
    // FXIME: Ts = Value or Struct
    template<typename ...Ts>
    struct Struct
    {
        using Type = ctll::list<Ts...>;
        static constexpr inline auto type = Type{};

        static constexpr std::size_t size() noexcept { return detail::size(type); }
    };
    */

    // FXIME: T = Value or Struct
    template<typename T, typename ArraySize, int N = 1>
    struct Array
    {
        using Type = T;
        using Size = ArraySize;
        static constexpr inline auto type = Type{};
        static constexpr inline auto n = N;

        static constexpr std::size_t stride() noexcept
        { return sizeOf(T{}); }
    };

    struct ArraySize : Value<SizeT> {};

    template<typename ArraySize, int N = 1>
    struct BitArray
    {
        //using Type = bool;
        using Size = ArraySize;
        //static constexpr inline auto type = Type{};
        static constexpr inline auto n = N;

        //static constexpr std::size_t stride() noexcept { return T::size(); }

        static constexpr auto byteIndex(SizeT index) noexcept
        {
            return index / 8;
        }

        static constexpr auto bitIndex(SizeT index) noexcept
        {
            return index % 8;
        }

        static constexpr auto bitMask(SizeT index) noexcept
        {
            return std::byte{1} << bitIndex(index);
        }


        static constexpr auto size(SizeT numBits) noexcept
        {
            auto n = byteIndex(numBits);
            if (bitIndex(numBits))
                ++n;
            return n;
        }

        static constexpr auto get(std::byte *data, SizeT index) noexcept
        {
            return (*(data + byteIndex(index)) & bitMask(index)) != std::byte{0};
        }

        static constexpr auto set(std::byte *data, SizeT index) noexcept
        {
            *(data + byteIndex(index)) |= bitMask(index);
        }

        static constexpr auto unset(std::byte *data, SizeT index) noexcept
        {
            *(data + byteIndex(index)) &= ~bitMask(index);
        }

    };

    struct BitValue : Value<bool> {}; // FIXME : add T to BitArray

    /*
    template<typename T> constexpr auto isValue(Value<T>) { return true; }
    constexpr auto isValue(...) { return false; }
    */

    /*
    template<typename ...Ts> constexpr auto isStruct(Struct<Ts...>) { return true; }
    constexpr auto isStruct(...) { return false; }
    */

    template<typename T, typename S, int N> constexpr auto isArray(Array<T, S, N>) { return true; }
    constexpr auto isArray(...) { return false; }

    template<typename S, int N> constexpr auto isBitArray(BitArray<S, N>) { return true; }
    constexpr auto isBitArray(...) { return false; }

} // namespace Kitimar::CTLayout
