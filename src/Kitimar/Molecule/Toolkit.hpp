#pragma once

namespace Kitimar::Toolkit {

    using Id = std::size_t;

    inline constexpr Id toolkitId(std::string_view name)
    {
        auto charValue = [] (char c) {
            if (c >= '0' && c <= '9')
                return c - '0';
            if (c >= 'A' && c <= 'Z')
                return c - 'A' + 10;
            if (c >= 'a' && c <= 'z')
                return c - 'a' + 10 + 26;
            if (c == '_')
                return 10 + 26 + 26;
            throw std::runtime_error("Toolkit name should only contain alphanumeric characters and '_'");
        };

        if (name.size() > 10)
            throw std::runtime_error("Toolkit name should be 10 characters or less");

        std::size_t h = 0;
        for (auto c : name) {
            h = (h << 6) + charValue(c);
        }
        return h;
    }

    template<Id ToolkitId>
    auto readSmiles(std::string_view smiles);

    template<Id ToolkitId, typename Mol>
    auto writeSmiles(const Mol &mol);

    template<Id ToolkitId>
    auto smilesMolSource(std::string_view filename);

    template<Id ToolkitId, typename Mol>
    auto match(std::string_view SMARTS, const Mol &mol);

    template<Id ToolkitId, typename Mol>
    auto map(std::string_view SMARTS, const Mol &mol);

    template<Id ToolkitId, typename Mol>
    auto count_unique(std::string_view SMARTS, const Mol &mol);

    template<Id ToolkitId, typename Mol>
    auto count_all(std::string_view SMARTS, const Mol &mol);

    template<Id ToolkitId, typename Mol>
    auto maps_unique(std::string_view SMARTS, const Mol &mol);

    template<Id ToolkitId, typename Mol>
    auto maps_all(std::string_view SMARTS, const Mol &mol);




} // namespace Kitimar::Toolkit
