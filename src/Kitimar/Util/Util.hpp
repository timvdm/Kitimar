#pragma once

#include <ctll/fixed_string.hpp>

namespace Kitimar::Util {



    template<auto N>
    auto toString(ctll::fixed_string<N> str)
    {
        return std::string{str.begin(), str.end()};
    }



} // namespace Kitimar::Util
