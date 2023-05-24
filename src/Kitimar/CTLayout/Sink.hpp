#pragma once

#include <Kitimar/CTLayout/Source.hpp>

#include <type_traits>
#include <concepts>
#include <string>

#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>

namespace Kitimar::CTLayout {

    template<typename Ptr = std::byte*>
    class PtrSink : public PtrSource<Ptr>
    {
        public:
            constexpr PtrSink(Ptr data = nullptr) noexcept : PtrSource<Ptr>{data}
            {
            }

            template<typename T>
            constexpr auto write(T value, std::size_t offset = 0) const noexcept
            {
                *reinterpret_cast<T*>(this->m_data + offset) = value;
            }

            constexpr auto operator+(std::size_t offset) const noexcept
            {
                return PtrSink<Ptr>(this->m_data + offset);
            }
    };

    using BytePtrSink = PtrSink<>;

    class StlVectorSink
    {
        public:
            constexpr StlVectorSink(std::vector<std::byte> &data, std::size_t offset = 0) : m_data{data}, m_offset(offset)
            {
            }

            constexpr auto& data() { return m_data; }
            constexpr const auto& data() const { return m_data; }

            constexpr operator bool() const noexcept
            {
                return true;
            }

            void resize(std::size_t size)
            {
                m_data.resize(size);
            }

            auto offset() const
            {
                return m_offset;
            }

            template<typename T>
            void write(T value, std::size_t offset = 0)
            {
                //std::cout << "StlVectorSink::write(value = " << value << ", offset = " << m_offset + offset << ", size = " << sizeof(T) << ")" << std::endl;
                if (m_data.size() < m_offset + offset + sizeof(T))
                    resize(m_offset + offset + sizeof(T));
                *reinterpret_cast<T*>(m_data.data() + m_offset + offset) = value;
            }

            template<typename T>
            auto read(std::size_t offset = 0) const
            {
                //std::cout << "StlVectorSink::read(offset = " << m_offset + offset << ", size = " << sizeof(T) << ")" << std::endl;
                assert(m_offset + offset + sizeof(T) <= m_data.size());
                return *reinterpret_cast<const T*>(m_data.data() + m_offset + offset);
            }

            constexpr StlVectorSink operator+(std::size_t offset) const noexcept
            {
                return StlVectorSink{m_data, m_offset + offset};
            }


        private:
            std::vector<std::byte> &m_data;
            std::size_t m_offset = 0;
    };

    class FileStreamSink
    {
        public:
            FileStreamSink(std::string_view path)
            {
                m_fs = std::make_shared<std::fstream>(path.data(), std::ios_base::binary | std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
                assert(*m_fs.get());
            }

            FileStreamSink(const std::shared_ptr<std::fstream> &fs, std::size_t offset = 0) : m_fs{fs}, m_offset{offset}
            {
            }

            operator bool() const noexcept
            {
                return m_fs.get();
            }

            auto& fs() { return m_fs; }
            const auto& fs() const { return m_fs; }

            void close()
            {
                m_fs->close();
            }

            template<typename T>
            auto read(std::size_t offset = 0) const
            {
                //std::cout << "read(offset = " << offset << ")" << std::endl;
                assert(*m_fs.get());
                m_fs->seekg(m_offset + offset);
                assert(*m_fs.get());
                T value;
                m_fs->read(reinterpret_cast<char*>(&value), sizeof(T));
                return value;
            }

            template<typename T>
            void write(T value, std::size_t offset = 0)
            {
                //std::cout << "write(value = " << value << ", offset = " << m_offset + offset << ")" << std::endl;
                assert(*m_fs.get());
                m_fs->seekp(m_offset + offset);
                assert(*m_fs.get());
                m_fs->write(reinterpret_cast<const char*>(&value), sizeof(T));
            }

            FileStreamSink operator+(std::size_t offset) const noexcept
            {
                return FileStreamSink{m_fs, m_offset + offset};
            }

        private:
            std::shared_ptr<std::fstream> m_fs;
            std::size_t m_offset = 0;
    };

} // namespace Kitimar::CTLayout
