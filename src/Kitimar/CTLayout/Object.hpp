#pragma once

#include <Kitimar/CTLayout/Layout.hpp>

namespace Kitimar::CTLayout {

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


} // namespace Kitimar::CTLayout
