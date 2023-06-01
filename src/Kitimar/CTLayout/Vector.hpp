#pragma once

#include <Kitimar/CTLayout/Source.hpp>

#include <limits>
#include <cassert>
#include <type_traits>
#include <cstdint> // uintX_t
#include <concepts>
#include <algorithm>
#include <vector>

namespace Kitimar::CTLayout {

    // FXIME: T = integral or floating_point
    template<typename T, typename LengthT = uint32_t>
    struct Vector
    {
        using Type = T;
        using Length = LengthT;
        using Index = LengthT;

        static constexpr inline auto type = Type{};
        static constexpr inline auto length = Length{};

        static constexpr inline auto invalidIndex = std::numeric_limits<Length>::max();
    };

    template<typename T>
    concept IsVector = std::derived_from<T, Vector<typename T::Type, typename T::Length>>;

    // isVector

    template<typename T, typename Length>
    constexpr bool isVector(Vector<T, Length>) noexcept { return true; }
    constexpr bool isVector(...) noexcept { return false; }    

    // isFixedSize

    template<typename T, typename Length>
    constexpr bool isFixedSize(Vector<T, Length>) noexcept { return false; }

    // contains

    template<typename VectorT, typename T>
    constexpr bool contains(VectorT, T) noexcept requires IsVector<VectorT>
    {
        if constexpr (std::is_same_v<VectorT, T>)
            return true;
        return contains(VectorT::type, T{});
    }

    // sizeOf

    template<typename Length>
    constexpr auto elementOffset(auto data, Length length, auto index)
    {
        if (!index)
            return sizeof(Length) + length * sizeof(std::size_t);
        auto off = sizeof(Length) + index * sizeof(std::size_t);
        return data.template read<std::size_t>(off);
    }

    template<typename VectorT, typename Source = PtrSource>
    constexpr std::size_t sizeOf(VectorT, Source data) noexcept requires IsVector<VectorT>
    {
        using Length = typename VectorT::Length;                                                                                                      

        if constexpr (isFixedSize(VectorT::type)) {
            auto length = data.template read<Length>();
            return sizeof(Length) + length * sizeOf(VectorT::type);
        } else {
            return data.template read<std::size_t>(sizeof(Length));
        }
    }    

    // offset

    template<typename VectorT, typename T, typename Source = PtrSource, typename Index = VectorT::Index, typename ...Indexes>
    constexpr std::size_t offset(VectorT, T, Source data = {}, Index index = {}, Indexes ...indexes) noexcept requires IsVector<VectorT>
    {
        using Length = typename VectorT::Length;

        if constexpr (std::is_same_v<VectorT, T>)
            return 0;
        assert(data);
        if constexpr (!contains(VectorT{}, T{}))
            return sizeOf(VectorT{}, data);
        else if constexpr (isFixedSize(VectorT::type))
            return sizeof(Length) + offset(VectorT::type, T{}) + index * sizeOf(VectorT::type);
        else {
            auto length = data.template read<Length>();
            auto off = elementOffset<Length>(data, length, index);
            return off + offset(VectorT::type, T{}, data + off, indexes...);
        }
    }

    // VectorObject

    template<typename VectorT, typename Source = PtrSource>
    class VectorObject
    {
        public:
            using Type = VectorT::Type;
            using Length = VectorT::Length;

            constexpr VectorObject(VectorT = {}, Source data = {}) noexcept : m_data(data)
            {
            }

            constexpr Length length() const noexcept
            {
                return m_data.template read<Length>();
            }

            constexpr auto at(std::integral auto index) const noexcept
            {
                auto off = offset(VectorT{}, Type{}, m_data, index);
                return toObject(Type{}, m_data + off);
            }

            constexpr auto range() const noexcept
            {
                static_assert(isValue(Type{}));
                auto off = offset(VectorT{}, Type{}, m_data, 0);
                return m_data.template range<typename Type::Type>(length(), off);
            }

        private:
            Source m_data = {};
    };

    // toObject

    template<typename VectorT, typename Source = PtrSource>
    constexpr auto toObject(VectorT, Source data) noexcept requires IsVector<VectorT>
    {
        return VectorObject{VectorT{}, data};
    }

    // vectorWriter

    template<typename VectorT, typename Sink, typename Finished = std::nullptr_t>
    class VectorWriter
    {
        public:
            using Type = VectorT::Type;
            using Length = VectorT::Length;

            VectorWriter(VectorT, Sink sink, Finished finished = nullptr) : m_sink(sink), m_finished{finished}
            {
                //fmt::println("VectorWriter::VectorWriter(offset={})", sink.offset());
            }

            constexpr Length length() const noexcept
            {
                assert(m_written.size());
                //fmt::println("VectorWriter::length()");
                return m_sink.template read<Length>();
            }

            constexpr void setLength(Length length) noexcept
            {
                //fmt::println("VectorWriter::setLength({})", length);
                m_sink.write(length);
                if (!length) {
                    if constexpr (!std::is_same_v<Finished, std::nullptr_t>)
                        m_finished();
                } else {
                    m_written.resize(length);
                    m_locked.resize(length);
                }
            }

            constexpr auto at(std::integral auto index) noexcept
            {
                //fmt::println("VectorWriter::at({})", index);
                assert(m_written.size());
                assert(!m_locked[index]);
                m_locked[index] = true;
                if (index)
                    assert(std::find(m_written.begin(), m_written.begin() + index, false) == m_written.begin() + index);

                auto off = offset(VectorT{}, Type{}, m_sink, index);                
                return toWriter(Type{}, m_sink + off, [this, off, index] () {
                    m_written[index] = true;
                    if constexpr (!isFixedSize(Type{})) {
                        //fmt::println("update({})", index);
                        auto off2 = elementOffset(m_sink, length(),  index) + sizeOf(Type{}, m_sink + off);
                        if (index + 1 < length())
                            m_sink.write(static_cast<std::size_t>(off2), sizeof(Length) + (index + 1) * sizeof(std::size_t));
                        else
                            m_sink.write(static_cast<std::size_t>(off2), sizeof(Length));

                    }
                    if constexpr (!std::is_same_v<Finished, std::nullptr_t>)
                        if (std::ranges::count(m_written, true) == m_written.size())
                            m_finished();
                });
            }

        private:
            Sink m_sink;
            std::vector<bool> m_written;
            std::vector<bool> m_locked;
            Finished m_finished;
    };

    // toWriter

    template<typename VectorT, typename Sink, typename Finished = std::nullptr_t>
    constexpr auto toWriter(VectorT, Sink sink, Finished finished = nullptr) noexcept requires IsVector<VectorT>
    {
        return VectorWriter{VectorT{}, sink, finished};
    }

} // namespace Kitimar::CTLayout
