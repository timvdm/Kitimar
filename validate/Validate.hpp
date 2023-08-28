#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Util/Util.hpp>
#include <Kitimar/RDKit/RDKit.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>

#include <catch2/catch_all.hpp>

using namespace Kitimar;

template<Toolkit::Id ToolkitId, ctll::fixed_string SMARTS>
void toolkit_validate(auto &mol)
{
    std::size_t ctseCount = ctse::count_unique<SMARTS>(mol);
    std::size_t toolkitCount = Toolkit::count_unique<ToolkitId>(Util::toString(SMARTS), mol);
    if (ctseCount != toolkitCount)
        std::cout << "FAIL: " << Util::toString(SMARTS) << '\n';
    CHECK(ctseCount == toolkitCount);
}

template<ctll::fixed_string SMARTS>
void validate(OpenBabel::OBMol &mol)
{
    toolkit_validate<Toolkit::openbabel, SMARTS>(mol);
}

template<ctll::fixed_string SMARTS>
void validate(RDKit::ROMol &mol)
{
    toolkit_validate<Toolkit::rdkit, SMARTS>(mol);
}
