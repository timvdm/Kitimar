#pragma once

#include <array>
#include <vector>
#include <ranges>
#include <cassert>

#ifdef KITIMAR_WITH_IOSTREAM
#include <iostream>

template<std::integral I, auto N>
std::ostream& operator<<(std::ostream &os, const std::array<I, N> &map)
{
    os << "[";
    for (auto i = 0; i < map.size(); ++i)
        os << " " << map[i];
    os << " ]";
    return os;
}

template<std::integral I>
std::ostream& operator<<(std::ostream &os, const std::vector<I> &v)
{
    os << "[ ";
    for (auto i : v)
        os << i << " ";
    os << "]";
    return os;
}

#endif // KITIMAR_WITH_IOSTREAM


namespace Kitimar::CTSmarts {


    template<std::integral Index, int N>
    using IsomorphismMap = std::array<Index, N>;

    template<std::integral Index, int N>
    using IsomorphismMaps = std::vector<IsomorphismMap<Index, N>>;

    template<std::integral Index, int NumQueryAtoms>
    class LookupMap
    {
        public:
            static constexpr inline auto invalidIndex = static_cast<Index>(-1);
            using Map = IsomorphismMap<Index, NumQueryAtoms>;

            constexpr LookupMap() noexcept
            {
                m_map.fill(invalidIndex);
            }

            constexpr const Map& map() const noexcept
            {
                return m_map;
            }

            constexpr Index operator()(int queryAtomIndex) const noexcept
            {
                return m_map[queryAtomIndex];
            }

            constexpr bool contains(int queryAtomIndex, auto atomIndex) const noexcept
            {
                return m_map[queryAtomIndex] == atomIndex;
            }

            constexpr bool containsAtom(auto atomIndex) const noexcept
            {
                return std::ranges::find(m_map, atomIndex) != m_map.end();
            }

            constexpr bool containsQueryAtom(int queryAtomIndex) const noexcept
            {
                return m_map[queryAtomIndex] != invalidIndex;
            }

            constexpr void reset(auto numAtoms) const noexcept {}

            constexpr void add(int queryAtomIndex, auto atomIndex) noexcept
            {
                m_map[queryAtomIndex] = atomIndex;
            }

            constexpr void remove(int queryAtomIndex, auto atomIndex) noexcept
            {
                m_map[queryAtomIndex] = invalidIndex;
            }

        protected:
            Map m_map;
    };

    template<std::integral Index, int NumAtoms>
    class InverseMap : public LookupMap<Index, NumAtoms>
    {
        public:
            constexpr void reset(auto numAtoms)
            {
                m_mapped.clear();
                m_mapped.resize(numAtoms);
            }

            constexpr bool containsAtom(auto atomIndex) const noexcept
            {
                assert(atomIndex < m_mapped.size());
                return m_mapped[atomIndex];
            }

            constexpr void add(int queryAtomIndex, auto atomIndex) noexcept
            {
                LookupMap<Index, NumAtoms>::add(queryAtomIndex, atomIndex);
                assert(atomIndex < m_mapped.size());
                assert(!m_mapped[atomIndex]);
                m_mapped[atomIndex] = true;
            }

            constexpr void remove(int queryAtomIndex, auto atomIndex) noexcept
            {
                LookupMap<Index, NumAtoms>::remove(queryAtomIndex, atomIndex);
                assert(atomIndex < m_mapped.size());
                assert(m_mapped[atomIndex]);
                m_mapped[atomIndex] = false;
            }

        private:
            std::vector<uint8_t> m_mapped;
    };

} // namespace Kitimar::CTSmarts
