#pragma once

namespace Kitimar::CTSmarts {

    struct UnconditionalFilter
    {
        static consteval bool enable([[maybe_unused]] auto smarts) noexcept
        {
            return true;
        }
    };

} // namespace Kitimar::CTSmarts
