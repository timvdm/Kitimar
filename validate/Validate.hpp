#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Util/Util.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <openbabel/parsmart.h>
#include <gtest/gtest.h>

using namespace Kitimar;

template<ctll::fixed_string SMARTS>
void validate_contains(OpenBabel::OBMol &mol)
{
    auto ctMatch = CTSmarts::contains<SMARTS>(mol);

    OpenBabel::OBSmartsPattern obsmarts;
    obsmarts.Init(Util::toString(SMARTS));
    auto obMatch = obsmarts.Match(mol);

    if (ctMatch != obMatch)
        std::cout << "FAIL: " << Util::toString(SMARTS) << std::endl;
    EXPECT_EQ(ctMatch, obMatch);
}
