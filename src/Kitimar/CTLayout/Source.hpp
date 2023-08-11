#pragma once

#include <mio/mmap.hpp>

#include <cstddef>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <vector>
#include <ranges>
#include <iostream>
#include <cassert>

namespace Kitimar::CTLayout {

    template<typename Ptr = const std::byte*>
    class PtrSourceBase
    {
        public:
            constexpr PtrSourceBase(Ptr data = nullptr) noexcept : m_data{data}
            {
            }

            operator bool() const noexcept
            {
                return m_data;
            }

            template<typename T>
            constexpr auto read(std::size_t offset = 0) const noexcept
            {
                return *reinterpret_cast<const T*>(m_data + offset);
            }

            template<typename T>
            constexpr auto range(auto length, std::size_t offset = 0) const noexcept
            {
                auto ptr = reinterpret_cast<const T*>(m_data + offset);
                return std::ranges::subrange(ptr, ptr + length);
            }

        protected:
            Ptr m_data;
    };

    class PtrSource : public PtrSourceBase<>
    {
        public:
            constexpr PtrSource(const std::byte *data = nullptr) noexcept : PtrSourceBase<>{data}
            {
            }

            constexpr auto operator+(std::size_t offset) const noexcept
            {
                return PtrSource(m_data + offset);
            }
    };


    class FileStreamSource
    {
        public:
            FileStreamSource() {}

            FileStreamSource(std::string_view path)
            {                
                m_ifs = std::make_shared<std::ifstream>(path.data(), std::ios_base::binary);
                assert(m_ifs);                
            }

            FileStreamSource(const std::shared_ptr<std::ifstream> &ifs, std::size_t offset = 0) : m_ifs{ifs}, m_offset{offset}
            {
            }

            operator bool() const noexcept
            {
                return m_ifs.get();
            }

            template<typename T>
            auto read(std::size_t offset = 0) const
            {
                assert(m_ifs);
                m_ifs->seekg(m_offset + offset);
                T value;
                m_ifs->read(reinterpret_cast<char*>(&value), sizeof(T));
                return value;
            }

            template<typename T>
            constexpr auto range(auto length, std::size_t offset = 0) const noexcept
            {
                std::vector<T> items;
                for (auto i = 0; i < length; ++i)
                    items.push_back(read<T>(offset + i * sizeof(T)));
                return items;
            }

            FileStreamSource operator+(std::size_t offset) const noexcept
            {
                return FileStreamSource{m_ifs, m_offset + offset};
            }                    

            auto& ifs()
            {
                return m_ifs;
            }

            const auto& ifs() const
            {
                return m_ifs;
            }

            void close()
            {
                m_ifs->close();
            }


        private:
            std::shared_ptr<std::ifstream> m_ifs;
            std::size_t m_offset = 0;
    };


    class InMemorySource
    {
        public:
            InMemorySource() {}

            InMemorySource(std::string_view filename)
            {
                std::ifstream ifs{filename.data(), std::ios_base::binary | std::ios_base::in};
                assert(ifs);
                auto size = std::filesystem::file_size(filename);
                m_data = std::make_shared<std::vector<std::byte>>(size);
                std::transform(std::istreambuf_iterator<char>(ifs),
                               std::istreambuf_iterator<char>(),
                               m_data->begin(),
                               [] (auto c) { return static_cast<std::byte>(c); });                
            }

            InMemorySource(const std::shared_ptr<std::vector<std::byte>> &data, std::size_t offset = 0) : m_data{data}, m_offset{offset}
            {
            }

            operator bool() const noexcept
            {
                return m_data.get();
            }

            template<typename T>
            auto read(std::size_t offset = 0) const
            {
                assert(m_data);                
                assert(m_offset + offset + sizeof(T) <= m_data->size());
                return *reinterpret_cast<T*>(m_data->data() + m_offset + offset);
            }

            template<typename T>
            auto range(auto length, std::size_t offset = 0) const noexcept
            {
                assert(m_data);
                assert(m_offset + offset + length * sizeof(T) <= m_data->size());
                auto ptr = reinterpret_cast<const T*>(m_data->data() + m_offset + offset);
                return std::ranges::subrange(ptr, ptr + length);
            }

            PtrSource toPtrSource() const noexcept
            {
                assert(m_data);
                return {m_data->data()};
            }

            auto operator+(std::size_t offset) const noexcept
            {
                return InMemorySource{m_data, m_offset + offset};
            }

        private:            
            std::shared_ptr<std::vector<std::byte>> m_data;
            std::size_t m_offset = 0;
    };

    using MioMemMapSource = mio::basic_mmap_source<std::byte>;

    class MemoryMappedSource
    {
        public:
            MemoryMappedSource() {}

            MemoryMappedSource(std::string_view filename)
            {
                m_source = std::make_shared<MioMemMapSource>(filename);
                assert(m_source->is_mapped());
            }

            MemoryMappedSource(const std::shared_ptr<MioMemMapSource> &source, std::size_t offset = 0) : m_source{source}, m_offset{offset}
            {
            }

            operator bool() const noexcept
            {
                return m_source->is_mapped();
            }

            template<typename T>
            auto read(std::size_t offset = 0) const
            {
                assert(m_source->is_mapped());
                assert(m_offset + offset + sizeof(T) <= m_source->size());
                return *reinterpret_cast<const T*>(m_source->begin() + m_offset + offset);
            }

            template<typename T>
            constexpr auto range(auto length, std::size_t offset = 0) const noexcept
            {
                assert(m_source->is_mapped());
                assert(m_offset + offset + length * sizeof(T) <= m_source->size());
                auto ptr = reinterpret_cast<const T*>(m_source->begin() + m_offset + offset);
                return std::ranges::subrange(ptr, ptr + length);
            }

            PtrSource toPtrSource() const noexcept
            {
                assert(m_source->is_mapped());
                return {m_source->begin()};
            }

            auto operator+(std::size_t offset) const noexcept
            {
                return MemoryMappedSource{m_source, m_offset + offset};
            }

        private:
            std::shared_ptr<MioMemMapSource> m_source;
            std::size_t m_offset = 0;
    };

} // namespace Kitimar::CTLayout
