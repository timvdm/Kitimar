#include <Kitimar/CTSmarts/CTSmarts.hpp>
#include <Kitimar/Util/Util.hpp>
#include <Kitimar/OpenBabel/OpenBabel.hpp>
#include <openbabel/parsmart.h>
#include <catch2/catch_all.hpp>

using namespace Kitimar;

template<ctll::fixed_string SMARTS>
void validate_contains(OpenBabel::OBMol &mol)
{
    auto ctMatch = CTSmarts::match<SMARTS>(mol);

    OpenBabel::OBSmartsPattern obsmarts;
    obsmarts.Init(Util::toString(SMARTS));
    auto obMatch = obsmarts.Match(mol);

    if (ctMatch != obMatch)
        std::cout << "FAIL: " << Util::toString(SMARTS) << '\n';
    CHECK(ctMatch == obMatch);
}
