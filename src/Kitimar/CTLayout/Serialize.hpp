#pragma once

#include <Kitimar/CTLayout/Object.hpp>

#include <mio/mmap.hpp>

#include <filesystem>
#include <fstream>
#include <vector>


namespace Kitimar::CTLayout {


    class FileStreamSink
    {
        public:
            FileStreamSink(std::string_view filename) : m_ofs(filename.data(), std::ios_base::binary | std::ios_base::out)
            {
                assert(m_ofs);
                // Placeholder for number of objects
                m_ofs.seekp(CTLayout::LayoutSize::size());
            }

            ~FileStreamSink()
            {
                close();
            }

            auto size() const noexcept
            {
                return m_size;
            }

            auto& ofs() noexcept
            {
                return m_ofs;
            }

            void close()
            {
                // Write size at start of file
                if (m_ofs) {
                    m_ofs.seekp(0);
                    m_ofs.write(reinterpret_cast<const char*>(&m_size), sizeof(m_size));
                }
                m_ofs.close();
            }

            auto write(const std::byte *data, auto size)
            {
                m_ofs.write(reinterpret_cast<const char*>(data), size);             
            }

            auto write(const std::vector<std::byte> &data)
            {
                m_ofs.write(reinterpret_cast<const char*>(data.data()), data.size());                
            }

            auto writeObject(const std::byte *data, auto size)
            {
                write(data, size);
                ++m_size;
            }

            auto writeObject(const std::vector<std::byte> &data)
            {
                write(data);
                ++m_size;
            }

        protected:
            std::ofstream m_ofs;
            LayoutSize::Type m_size = 0;
    };

    template<typename Layout>
    class FileStreamSource
    {
        public:
            FileStreamSource(std::string_view filename) : m_ifs{filename.data(), std::ios_base::binary | std::ios_base::in}
            {
                assert(m_ifs);
                // Read the number of objects
                m_ifs.read(reinterpret_cast<char*>(&m_size), LayoutSize::size());
            }

            auto size() const noexcept
            {
                return m_size;
            }

            auto index() const noexcept
            {
                return m_index;
            }

            auto& ifs() noexcept
            {
                return m_ifs;
            }

            Object<Layout> read()
            {
//                return read(m_buffer);
                if (m_index + 1 > m_size || !m_ifs)
                    return {};

                // Read size of layout
                m_buffer.resize(LayoutSize::size());
                m_ifs.read(reinterpret_cast<char*>(m_buffer.data()), m_buffer.size());
                if (!m_ifs)
                    return {};

                // Resize buffer
                auto size = *reinterpret_cast<LayoutSize::Type*>(m_buffer.data());
                m_buffer.resize(size);

                // Read remaining layout data
                m_ifs.read(reinterpret_cast<char*>(m_buffer.data() + LayoutSize::size()), size - LayoutSize::size());
                if (!m_ifs)
                    return {};

                ++m_index;
                return Object<Layout>{m_buffer.data()};
            }

            auto objects()
            {
                return std::views::iota(static_cast<decltype(size())>(0), size()) |
                        std::views::transform([this] (auto i) {
                            return read();
                        });
            }

            auto enumerate()
            {
                return std::views::iota(static_cast<decltype(size())>(0), size()) |
                        std::views::transform([this] (auto i) {
                            auto [ok, mol] = read();
                            return std::make_tuple(ok, i, mol);
                        });
            }

        private:
            std::ifstream m_ifs;
            std::vector<std::byte> m_buffer;
            LayoutSize::Type m_size = 0;
            LayoutSize::Type m_index = 0;
    };

    template<typename Layout>
    class InMemorySource
    {
        public:
            InMemorySource(std::string_view filename)
            {                
                std::ifstream ifs{filename.data(), std::ios_base::binary | std::ios_base::in};
                assert(ifs);
                auto size = std::filesystem::file_size(filename);
                m_data.resize(size);
                std::transform(std::istreambuf_iterator<char>(ifs),
                               std::istreambuf_iterator<char>(),
                               m_data.begin(),
                               [] (auto c) { return static_cast<std::byte>(c); });
                // Read the number of objects
                m_size = *reinterpret_cast<const LayoutSize::Type*>(m_data.data());
            }

            auto size() const noexcept
            {
                return m_size;
            }

            auto index() const noexcept
            {
                return m_index;
            }

            auto offset() const noexcept
            {
                return m_offset;
            }

            auto data() noexcept
            {
                return m_data.data();
            }

            auto read()
            {
                using Ptr = decltype(m_data.data());
                auto data = m_data.data() + m_offset;
                auto size = *reinterpret_cast<const LayoutSize::Type*>(data);
                ++m_index;
                m_offset += size;

                return Object<Layout, Ptr>{data};
            }

            auto objects()
            {
                return std::views::iota(static_cast<decltype(size())>(0), size()) |
                        std::views::transform([this] (auto i) {
                            return read();
                        });
            }

            auto enumerate()
            {
                return std::views::iota(static_cast<decltype(size())>(0), size()) |
                        std::views::transform([this] (auto i) {
                            return std::make_tuple(i, read());
                        });
            }

        private:
            std::vector<std::byte> m_data;
            LayoutSize::Type m_size = 0;
            LayoutSize::Type m_index = 0;
            LayoutSize::Type m_offset = LayoutSize::size();

    };

    using MioMemMapSource = mio::basic_mmap_source<std::byte>;

    template<typename Layout>
    class MemoryMappedSource
    {
        public:
            MemoryMappedSource(std::string_view filename) : m_source{filename}
            {                
                assert(m_source.is_mapped());
                // Read the number of objects
                m_size = *reinterpret_cast<const LayoutSize::Type*>(m_source.begin());
            }

            auto size() const noexcept
            {
                return m_size;
            }

            auto index() const noexcept
            {
                return m_index;
            }

            auto offset() const noexcept
            {
                return m_offset;
            }

            auto& source() noexcept
            {
                return m_source;
            }

            auto read()
            {
                using Ptr = decltype(m_source.begin());
                auto data = m_source.begin() + m_offset;
                auto size = *reinterpret_cast<const LayoutSize::Type*>(data);
                ++m_index;
                m_offset += size;
                return Object<Layout, Ptr>{data};
            }

            auto objects()
            {
                return std::views::iota(static_cast<decltype(size())>(0), size()) |
                        std::views::transform([this] (auto i) {
                            return read();
                        });
            }

            auto enumerate()
            {
                return std::views::iota(static_cast<decltype(size())>(0), size()) |
                        std::views::transform([this] (auto i) {
                            return std::make_tuple(i, read());
                        });
            }

        private:
            MioMemMapSource m_source;
            LayoutSize::Type m_size = 0;
            LayoutSize::Type m_index = 0;
            LayoutSize::Type m_offset = LayoutSize::size();
    };








} // namespace Kitimar::CTLayout
