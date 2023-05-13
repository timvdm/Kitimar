#pragma once

#include <ctll/list.hpp>

#include <cstdint>
#include <utility>
#include <array>
#include <algorithm>
#include <ranges>
#include <iostream>

#include <cassert>

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


namespace Kitimar {

    using SizeT = uint32_t;

    namespace detail {

        template<typename T, typename ...Ts>
        constexpr std::size_t size(ctll::list<T, Ts...>) noexcept { return T::size() + size(ctll::list<Ts...>{}); }
        constexpr std::size_t size(ctll::list<>) noexcept { return 0; }

    }


    // FXIME: T = integral or floating_point
    template<typename T>
    struct Value
    {
        using Type = T;
        static constexpr inline auto type = Type{};

        static constexpr std::size_t size() noexcept { return sizeof(Type); }
    };

    // FXIME: Ts = Value or Struct
    template<typename ...Ts>
    struct Struct
    {
        using Type = ctll::list<Ts...>;
        static constexpr inline auto type = Type{};

        static constexpr std::size_t size() noexcept { return detail::size(type); }
    };

    // FXIME: T = Value or Struct
    template<typename T, typename ArraySize, int N = 1>
    struct Array
    {
        using Type = T;
        using Size = ArraySize;
        static constexpr inline auto type = Type{};
        static constexpr inline auto n = N;

        static constexpr std::size_t stride() noexcept { return T::size(); }
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




    template<typename T> constexpr auto isValue(Value<T>) { return true; }
    constexpr auto isValue(...) { return false; }

    template<typename ...Ts> constexpr auto isStruct(Struct<Ts...>) { return true; }
    constexpr auto isStruct(...) { return false; }

    template<typename T, typename S, int N> constexpr auto isArray(Array<T, S, N>) { return true; }
    constexpr auto isArray(...) { return false; }

    template<typename S, int N> constexpr auto isBitArray(BitArray<S, N>) { return true; }
    constexpr auto isBitArray(...) { return false; }


    namespace detail {

        struct EndValue : Value<void*> {};

        template<typename T, typename Container>
        constexpr std::size_t find(T, Container) noexcept;

        template<typename T>
        constexpr std::size_t findStruct(T, ctll::empty_list) noexcept { return 0; }

        template<typename T, typename U, typename ...Us>
        constexpr std::size_t findStruct(T, ctll::list<U, Us...>) noexcept
        {
            auto pos1 = find(T{}, U{});
            if (pos1 < U::size())
                return pos1;

            return pos1 + findStruct(T{}, ctll::list<Us...>{});
        }



        template<typename T, typename Container>
        constexpr std::size_t find(T, Container) noexcept
        {
            if constexpr (std::is_same_v<T, Container>)
                return 0;
            else if constexpr (isStruct(Container{}))
                return findStruct(T{}, Container::type);
            else if constexpr (isArray(Container{}))
                return find(T{}, Container::type);
            else return Container::size();
        }

        constexpr auto findArrays(ctll::empty_list)
        {
            return ctll::empty_list{};
        }

        template<typename T, typename ...Ts>
        constexpr auto findArrays(ctll::list<T, Ts...>)
        {
            if constexpr (isArray(T{}) || isBitArray(T{}))
                return ctll::push_front(T{}, findArrays(ctll::list<Ts...>{}));
            else
                return findArrays(ctll::list<Ts...>{});
        }

        constexpr auto findArraySizes(ctll::empty_list)
        {
            return ctll::empty_list{};
        }

        template<typename T, typename ...Ts>
        constexpr auto findArraySizes(ctll::list<T, Ts...>)
        {
            if constexpr (isArray(T{}) || isBitArray(T{})) {
                auto arrays = findArraySizes(ctll::list<Ts...>{});
                if constexpr (ctll::exists_in(typename T::Size{}, arrays))
                    return arrays;
                else
                    return ctll::push_front(typename T::Size{}, arrays);
            } else
                return findArraySizes(ctll::list<Ts...>{});
        }


        template<typename T, typename U, typename ...Us>
        constexpr auto indexOf(T, ctll::list<U, Us...>) noexcept
        {
            if constexpr (std::is_same_v<T, U>)
                return 0;
            else
                return 1 + indexOf(T{}, ctll::list<Us...>{});
        }

    } // namespace detail

    struct LayoutSize : Value<std::size_t> {};

    template<typename ...Ts>
    struct Layout
    {
        using Type = ctll::list<Ts...>; // FIXME: unique + no array sizes
        static constexpr inline auto type = Type{};

        static constexpr inline auto arrays = detail::findArrays(type);
        static constexpr inline SizeT numArrays = ctll::size(arrays);

        static constexpr inline auto arraySizes = detail::findArraySizes(type);
        static constexpr inline auto numArraySizes = ctll::size(arraySizes);


        static constexpr inline auto header = ctll::push_front(LayoutSize{}, arraySizes);
        static constexpr inline auto headerSize = LayoutSize::size() + ctll::size(arraySizes) * sizeof(SizeT);
        using Header = std::remove_const_t<decltype(header)>;

        struct Loc
        {
            using ArraySizes = std::array<SizeT, numArraySizes>;

            std::size_t fixed = 0;
            ArraySizes arraySizes = {};
            ArraySizes bitArraySizes = {};
            SizeT stride = 0;

            constexpr auto offset(const SizeT *s) const noexcept
            {
                auto o = fixed;
                for (auto i = 0; i < arraySizes.size(); ++i) {
                    o += arraySizes[i] * s[i];
                    o += BitArray<ArraySize>::size(bitArraySizes[i] * s[i]);
                }
                return o;
            }

            constexpr auto offset(const SizeT *s, SizeT index) const noexcept
            {
                return offset(s) + index * stride;
            }

            friend constexpr auto operator+(const Loc &a, const Loc &b) noexcept
            {
                auto newArraySizes = a.arraySizes;
                for (auto i = 0; i < newArraySizes.size(); ++i) // FIXME: move to function...
                    newArraySizes[i] += b.arraySizes[i];
                auto newBitArraySizes = a.bitArraySizes;
                for (auto i = 0; i < newBitArraySizes.size(); ++i)
                    newBitArraySizes[i] += b.bitArraySizes[i];
                return Loc{ a.fixed + b.fixed, {newArraySizes}, {newBitArraySizes}, a.stride + b.stride };
            }

            friend constexpr std::ostream& operator<<(std::ostream &os, const Loc &l) noexcept
            {
                os << "Loc( offset = " << l.fixed << ", arraySizes = [ ";
                for (auto s : l.arraySizes)
                    os << s << " ";
                os << "], bitArraySizes = [ ";
                for (auto s : l.bitArraySizes)
                    os << s << " ";
                os << "], stride = " << l.stride << " )";
                return os;
            }
        };



        template<typename T>
        static constexpr auto _find(T, ctll::empty_list) noexcept
        {
            return Loc{};
        }

        template<typename T, typename U, typename ...Us>
        static constexpr auto _find(T, ctll::list<U, Us...>) noexcept
        {
            if constexpr (std::is_same_v<T, U>)
                return Loc{};
            else if constexpr (isArray(U{})) {
                auto size = detail::find(T{}, U{});
                if (size != U::stride())
                    return Loc{ size, {}, {}, U::stride() };

                auto i = detail::indexOf(typename U::Size{}, arraySizes);
                auto l = Loc{};
                l.arraySizes[i] = U::n * U::stride();
                return Loc{ 0UL, l.arraySizes} + _find(T{}, ctll::list<Us...>{}); // FIXME l + ...
            } else if constexpr (isBitArray(U{})) {
                auto i = detail::indexOf(typename U::Size{}, arraySizes);
                auto l = Loc{};
                l.bitArraySizes[i] = U::n;
                return Loc{ 0UL, {}, l.bitArraySizes } + _find(T{}, ctll::list<Us...>{}); // FIXME l + ...
            } else {
                auto size = detail::find(T{}, U{});
                if (size != U::size())
                    return Loc{ size };

                return Loc{ U::size() } + _find(T{}, ctll::list<Us...>{});
            }
        }

        template<typename T>
        static constexpr auto find(T) noexcept
        {
            return _find(T{}, type);
        }




        static constexpr auto size(const SizeT *sizes) noexcept
        {
            auto addr = find(detail::EndValue{});
            return headerSize + addr.offset(sizes);
        }

        static constexpr auto size(const Loc::ArraySizes &sizes) noexcept
        {
            return size(sizes.data());
        }

        template<typename T>
        static constexpr auto offset(const SizeT *sizes) noexcept
        {
            constexpr auto addr = find(T{});
            static_assert(!addr.stride, "T is an Array element, use offset(sizes, index)");
            return headerSize + addr.offset(sizes);
        }

        template<typename T>
        static constexpr auto offset(const Loc::ArraySizes &sizes) noexcept
        {
            return offset<T>(sizes.data());
        }

        template<typename T>
        static constexpr auto offset(const SizeT *sizes, SizeT index) noexcept
        {
            constexpr auto addr = find(T{});
            //static_assert(addr != size(sizes), "T is not an Array element, use offset(sizes)");
            static_assert(addr.stride, "T is not an Array element, use offset(sizes)");
            return headerSize + addr.offset(sizes, index);
        }

        template<typename T>
        static constexpr auto offset(const Loc::ArraySizes &sizes, SizeT index) noexcept
        {
            return offset<T>(sizes.data(), index);
        }

    };



    template<typename ObjectLayout, typename Ptr = const std::byte*> // FIXME: rename to It, random access iterator concept....
    class Object
    {
        public:
            using Layout = ObjectLayout;

            constexpr Object(Ptr data = nullptr) noexcept : m_data(data)
            {
            }

            constexpr auto data() const noexcept
            {
                return m_data;
            }

            constexpr auto setData(Ptr data) noexcept
            {
                m_data = data;
            }

            constexpr auto size() const noexcept
            {
                return Layout::size(sizes());
            }

            constexpr auto setSize(std::size_t value) const noexcept
            {
                auto &ref = *reinterpret_cast<std::size_t*>(m_data);
                ref = value;
            }

            template<typename A>
            constexpr auto size() const noexcept
            {
                static_assert(isArray(A{}) || isBitArray(A{}));
                auto index = detail::indexOf(typename A::Size{}, Layout::arraySizes);
                return *reinterpret_cast<const SizeT*>(m_data + LayoutSize::size() + index * sizeof(SizeT)) * A::n;
            }

            template<typename S>
            constexpr auto setSize(S, SizeT value) noexcept
            {
                static_assert(isValue(S{}));
                constexpr auto index = detail::indexOf(S{}, Layout::arraySizes);
                static_assert(index < ctll::size(Layout::arraySizes));
                auto &ref = *reinterpret_cast<SizeT*>(m_data + LayoutSize::size() + index * sizeof(SizeT));
                ref = value;
            }

            template<typename T>
            constexpr auto offset() const noexcept
            {
                return Layout::template offset<T>(sizes());
            }

            template<typename T>
            constexpr auto offset(SizeT index) const noexcept
            {
                return Layout::template offset<T>(sizes(), index);
            }

            template<typename V>
            constexpr auto value() const noexcept
            {
                static_assert(isValue(V{}));
                return *reinterpret_cast<const V::Type*>(m_data + offset<V>());
            }

            template<typename V>
            constexpr auto value(SizeT index) const noexcept
            {
                static_assert(isValue(V{}));
                return *reinterpret_cast<const V::Type*>(m_data + offset<V>(index));
            }

            template<typename V>
            constexpr auto setValue(V::Type value) noexcept
            {
                static_assert(isValue(V{}));
                auto &ref = *reinterpret_cast<V::Type*>(m_data + offset<V>());
                ref = value;
            }

            template<typename V>
            constexpr auto setValue(SizeT index, V::Type value) noexcept
            {
                static_assert(isValue(V{}));
                auto &ref = *reinterpret_cast<V::Type*>(m_data + offset<V>(index));
                ref = value;
            }

            template<typename V>
            constexpr auto setValue(V, SizeT index, V::Type value) noexcept
            {
                static_assert(isValue(V{}));
                auto &ref = *reinterpret_cast<V::Type*>(m_data + offset<V>(index));
                ref = value;
            }

            template<typename BA>
            constexpr auto bit(BA, SizeT index) noexcept
            {
                static_assert(isBitArray(BA{}));
                auto data = m_data + offset<BA>();
                return BA::get(data, index);
            }

            template<typename BA>
            constexpr auto setBit(BA, SizeT index) noexcept
            {
                static_assert(isBitArray(BA{}));
                auto data = m_data + offset<BA>();
                return BA::set(data, index);
            }

            template<typename BA>
            constexpr auto unsetBit(BA, SizeT index) noexcept
            {
                static_assert(isBitArray(BA{}));
                auto data = m_data + offset<BA>();
                return BA::unset(data, index);
            }

        private:
            constexpr auto sizes() const noexcept
            {
                return reinterpret_cast<const SizeT*>(m_data + LayoutSize::size());
            }

            Ptr m_data = nullptr;
    };

    template<typename ObjectLayout>
    using EditableObject = Object<ObjectLayout, std::byte*>;


} // namespace Kitimar
