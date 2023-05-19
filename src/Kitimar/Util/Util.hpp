#pragma once

#include <ctll/fixed_string.hpp>

#include <any>
#include <string>
#include <sstream>
#include <typeindex>
#include <functional>
#include <unordered_map>
#include <iostream>

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

    inline std::string toString(const std::any &value)
    {
        if (detail::anyToStringMap.contains(std::type_index(value.type())))
            return detail::anyToStringMap[std::type_index(value.type())](value);
        return "???";
    }

    template<auto N>
    std::string toString(ctll::fixed_string<N> str)
    {
        return std::string{str.begin(), str.end()};
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



} // namespace Kitimar::Util
