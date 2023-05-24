#pragma once

#include <mio/mmap.hpp>

#include <cstddef>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <vector>
//#include <algorithm>
#include <ranges>
#include <iostream>
#include <cassert>

namespace Kitimar::CTLayout {

    template<typename Ptr = const std::byte*>
    class PtrSource
    {
        public:
            constexpr PtrSource(Ptr data = nullptr) noexcept : m_data{data}
            {
            }

            constexpr operator Ptr() const noexcept
            {
                return m_data;
            }

            template<typename T>
            constexpr auto read(std::size_t offset = 0) const noexcept
            {
                return *reinterpret_cast<const T*>(m_data + offset);
            }

            constexpr auto operator+(std::size_t offset) const noexcept
            {
                return PtrSource<Ptr>(m_data + offset);
            }

        protected:
            Ptr m_data;
    };

    using BytePtrSource = PtrSource<>;

    class FileStreamSource
    {
        public:
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

            auto& ifs() { return m_ifs; }
            const auto& ifs() const { return m_ifs; }

            void close()
            {
                m_ifs->close();
            }

            template<typename T>
            auto read(std::size_t offset = 0) const
            {
                assert(m_ifs);
                m_ifs->seekg(offset);
                T value;
                m_ifs->read(reinterpret_cast<char*>(&value), sizeof(T));
                return value;
            }

            FileStreamSource operator+(std::size_t offset) const noexcept
            {
                return FileStreamSource{m_ifs, m_offset + offset};
            }



        private:
            std::shared_ptr<std::ifstream> m_ifs;
            std::size_t m_offset = 0;
    };


    class InMemorySource
    {
        public:
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
                //std::cout << "StlVectorSink::read(offset = " << m_offset + offset << ", size = " << sizeof(T) << ")" << std::endl;
                assert(m_offset + offset + sizeof(T) <= m_data->size());
                return *reinterpret_cast<T*>(m_data->data() + m_offset + offset);
            }

            operator const std::byte*() const noexcept
            {
                return m_data->data() + m_offset;
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
            MemoryMappedSource(std::string_view filename)
            {
                std::cout << filename << std::endl;
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
                //std::cout << "StlVectorSink::read(offset = " << m_offset + offset << ", size = " << sizeof(T) << ")" << std::endl;
                assert(m_offset + offset + sizeof(T) <= m_source->size());
                return *reinterpret_cast<const T*>(m_source->begin() + m_offset + offset);
            }

            operator const std::byte*() const noexcept
            {
                return m_source->begin() + m_offset;
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
