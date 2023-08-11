#pragma once

namespace Kitimar::CTSmarts {

    struct UnconditionalFilter
    {
        static consteval bool enable(auto smarts) noexcept
        {
            return true;
        }
    };

} // namespace Kitimar::CTSmarts
