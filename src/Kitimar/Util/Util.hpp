#pragma once

#include <any>
#include <string>
#include <sstream>
#include <typeindex>
#include <functional>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cassert>
#include <filesystem>
#include <algorithm>

#ifndef _MSC_VER
#include <cxxabi.h> // FIXME
#endif

namespace Kitimar::Util {

    namespace detail {

        using AnyToString = std::function<std::string(const std::any&)>;
        using AnyToStringMap = std::unordered_map<std::type_index, AnyToString>;


        template<typename T>
        constexpr auto anyIntegralToString() noexcept
        {
            return [] (const std::any &value) {
                return std::to_string(std::any_cast<T>(value));
            };
        }

        static inline AnyToStringMap anyToStringMap = {
            { std::type_index(typeid(void)), [] (const std::any &value) { return std::string{}; } },
            { std::type_index(typeid(int)), anyIntegralToString<int>() },
            { std::type_index(typeid(unsigned long)), anyIntegralToString<unsigned long>() }
        };

    } // namespace detail

    template<typename T>
    std::string typeName(T, bool qualified = false)
    {
        #ifndef _MSC_VER
        std::string_view name = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, nullptr);
        #else
        std::string_view name = typeid(T).name();
        #endif
        if (qualified)
            return name.data();
        auto pos = name.rfind("::");
        if (pos == std::string_view::npos)
            return name.data();
        return name.substr(pos + 2).data();
    }

    template<typename Container>
    auto indexOf(const Container &c, auto value)
    {
        return std::ranges::find(c, value) - c.begin();
    }


    inline std::string toString(const std::any &value)
    {
        if (detail::anyToStringMap.contains(std::type_index(value.type())))
            return detail::anyToStringMap[std::type_index(value.type())](value);
        return "???";
    }

    inline std::string pad(const std::string &str, std::size_t width = 30)
    {
        if (str.size() >= width)
            return str;
        return str + std::string(width - str.size(), ' ');
    }

    std::string pad(std::integral auto n, std::size_t width = 10)
    {
        std::stringstream ss;
        ss << n;
        ss << std::string(width - ss.str().size(), ' ');
        return ss.str();
    }

    inline std::vector<std::byte> readFileData(std::string_view path)
    {
        std::ifstream ifs{path.data(), std::ios_base::binary | std::ios_base::in};
        assert(ifs);
        auto size = std::filesystem::file_size(path);
        std::vector<std::byte> data;
        data.reserve(size);
        std::istreambuf_iterator<char> it{ifs}, end;
        for (; it != end; ++it)
            data.push_back(static_cast<std::byte>(*it));
        return data;
    }

    inline void writeFileData(std::string_view path, const std::vector<std::byte> &data)
    {
        std::ofstream ofs(path.data(), std::ios_base::binary | std::ios_base::out);
        assert(ofs);
        ofs.write(reinterpret_cast<const char*>(data.data()), data.size());
        ofs.close();
    }



} // namespace Kitimar::Util
