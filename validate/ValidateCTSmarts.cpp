#include "Validate.hpp"
#include "../test/TestData.hpp"
#include <catch2/catch_all.hpp>

using namespace Kitimar;

//
// Basic
//

#ifdef KITIMAR_WITH_VALIDATION_BASIC

template<Molecule::Molecule Mol>
void validate_basic(Mol &mol);

template<Toolkit::Id ToolkitId>
void toolkit_validate_basic()
{
    auto source = Toolkit::smilesMolSource<ToolkitId>(chembl_smi_filename("100K"));
    auto i = 0;
    for (auto mol : source.molecules()) {
        std::cout << "Molecule: " << ++i << " -- " << Toolkit::writeSmiles<ToolkitId>(*mol) << '\n';
        validate_basic(*mol);
    }
}

TEST_CASE("rdkit_basic")
{
    toolkit_validate_basic<Toolkit::rdkit>();
}


TEST_CASE("openbabel_basic")
{
    toolkit_validate_basic<Toolkit::openbabel>();
}

#endif // KITIMAR_WITH_VALIDATION_BASIC

//
// Substructure Query Collection
//

#ifdef KITIMAR_WITH_VALIDATION_SQC

template<Molecule::Molecule Mol>
void validate_parts(Mol &mol);

template<Toolkit::Id ToolkitId>
void toolkit_validate_parts()
{
    auto source = Toolkit::smilesMolSource<ToolkitId>(chembl_smi_filename());
    auto i = 0;
    for (auto mol : source.molecules()) {
        std::cout << "Molecule: " << ++i << " -- " << Toolkit::writeSmiles<ToolkitId>(*mol) << '\n';
        validate_parts(*mol);
    }
}

TEST_CASE("rdkit_sqc")
{
    toolkit_validate_parts<Toolkit::rdkit>();
}

TEST_CASE("openbabel_sqc")
{
    toolkit_validate_parts<Toolkit::openbabel>();
}

#endif // KITIMAR_WITH_VALIDATION_SQC

template<Toolkit::Id ToolkitId, ctll::fixed_string SMARTS>
void validate(const char *SMILES)
{
    std::cout << Toolkit::name<ToolkitId>() << ": " << SMILES << std::endl;
    auto mol = Toolkit::readSmiles<ToolkitId>(SMILES);
    validate<SMARTS>(*mol);
}

TEST_CASE("not_hydrogen")
{
    validate<Toolkit::openbabel, "[!H]">("C");
    validate<Toolkit::openbabel, "[!H]">("[H]"); // OpenBabel counts 1
    validate<Toolkit::openbabel, "[!H]">("c1ccccc1"); // OpenBabel counts 0

    validate<Toolkit::rdkit, "[!H]">("C");
    validate<Toolkit::rdkit, "[!H]">("[H]"); // RDKit counts 1
    validate<Toolkit::rdkit, "[!H]">("c1ccccc1"); // RDKit counts 6
}

TEST_CASE("rdkit_sqc_fails")
{
    //validate<Toolkit::rdkit, "[+,++,+++]">("[Cl+3]"); // “[+++]” not supported? (not standard)
    //validate<Toolkit::rdkit, "[-,--,---]">("[Al-3]"); // “[---]” not supported? (not standard)
    //validate<Toolkit::rdkit, "[+H]">("C[NH+](C)C"); // RDkit counts 1


    // aromatic triple bond (FIXED)
    validate<Toolkit::rdkit, "*1****1">("c1nc#cn1C");
    validate<Toolkit::rdkit, "*1*****1">("c1ccccc1");
    validate<Toolkit::rdkit, "a1aaaa1">("c1nc#cn1C");
    validate<Toolkit::rdkit, "[#6]#[#6]">("c1nc#cn1C");
    validate<Toolkit::rdkit, "*-*-[R]~*:[a]">("CCc1nc#cn1CC1CCC(=O)O1");
    validate<Toolkit::rdkit, "*:[a]">("CCc1nc#cn1CC1CCC(=O)O1");
    validate<Toolkit::rdkit, "*:*">("CCc1nc#cn1CC1CCC(=O)O1");
    validate<Toolkit::rdkit, "*-*-[R]~*~*:[a]">("CCc1nc#cn1CC1CC(c2ccccc2)(c2ccccc2)C(=O)O1");

    // ???
    validate<Toolkit::rdkit, "*!:[a]:*:*:[a]!:*">("");
    validate<Toolkit::rdkit, "[!#1]!:*:*!:[!#1]">("");
    validate<Toolkit::rdkit, "[#6]:[#6$(a1aaaa1)]">("");
    validate<Toolkit::rdkit, "[#6]:[#7$(a1aaaa1)]">("");
    validate<Toolkit::rdkit, "[R](-*(-*))~*~*~*~[a]">("");
    validate<Toolkit::rdkit, "[#9,#17,#35,#53]~*(~*)~*">("");
    validate<Toolkit::rdkit, "[F,Cl,Br,I]~*(~[!#1])~[!#1]">("");
    validate<Toolkit::rdkit, "[!#1]~*(~[!#1])(~[!#1])~[!#1]">("");
    validate<Toolkit::rdkit, "[!#1]~[!#6;!#1](~[!#1])~[!#1]">("");
    validate<Toolkit::rdkit, "[$([cX3](:*):*),$([cX2+](:*):*)]">("");
    validate<Toolkit::rdkit, "[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]">("");
    validate<Toolkit::rdkit, "[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]">("");
    validate<Toolkit::rdkit, "[$([cX3](:*):*),$([cX2+](:*):*),$([CX3]=*),$([CX2+]=*)]">("");
    validate<Toolkit::rdkit, "[Be,B,Al,Ti,Cr,Mn,Fe,Co,Ni,Cu,Pd,Ag,Sn,Pt,Au,Hg,Pb,Bi,As,Sb,Gd,Te]">("");

    // !@
    //validate<Toolkit::rdkit, "*!@*">("CCCCCCCCCC(C(=O)NCCc1ccc(OP(=S)(Oc2ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc2)N(C)/N=C/c2ccc(OP3(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)=NP(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)=NP(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)=N3)cc2)cc1)P(=O)(O)O.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO");
    //validate<Toolkit::rdkit, "*(!@*)(!@*)">("CCCCCCCCCC(C(=O)NCCc1ccc(OP(=S)(Oc2ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc2)N(C)/N=C/c2ccc(OP3(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)=NP(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)=NP(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)(Oc4ccc(/C=N/N(C)P(=S)(Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)Oc5ccc(CCNC(=O)C(CCCCCCCCC)P(=O)(O)O)cc5)cc4)=N3)cc2)cc1)P(=O)(O)O.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO.CCCCCCCCCCCCCCCCNOC(CO)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)CO");

    //validate<Toolkit::rdkit, "[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]~[AR0]">("CC[C@H](C)[C@H](NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)NC(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](N)CCCNC(=N)N)[C@@H](C)CC)[C@@H](C)CC)C(=O)NCC(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)N[C@@H](CS)C(=O)O");
    //validate<Toolkit::rdkit, "[!#6]~*~*~[R]">("C=C1NC(=O)C(C)=CN1[C@H]1C[C@H](OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(=O)[nH]c(N)nc43)[C@H](O)[C@@H]2OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2OP(=O)(O)OC[C@H]2O[C@@H](n3ccc(N)nc3=O)[C@H](O)[C@@H]2OP(=O)(O)OC[C@H]2O[C@@H](n3cc(C)c(=O)[nH]c3=O)C[C@@H]2OP(=O)(O)OC[C@]23CNO[C@H]([C@@H]2O)[C@H](n2cc(C)c(=O)[nH]c2=O)O3)[C@@H](COP(=O)(O)O[C@@H]2[C@@H](COP(=O)(O)O[C@@H]3[C@@H](COP(=O)(O)O[C@@H]4[C@@H](COP(=O)(O)O[C@@H]5[C@@H](COP(=O)(O)O[C@H]6C[C@H](n7cc(C)c(=O)[nH]c7=O)O[C@@H]6COP(=O)(O)O[C@H]6C[C@H](n7cc(C)c(=O)[nH]c7=O)O[C@@H]6COP(=O)(O)O[C@@H]6[C@@H](COP(=O)(O)O[C@@H]7[C@@H](COP(=O)(O)O[C@@H]8[C@@H](COP(=O)(O)O[C@@H]9[C@@H](COP(=O)(O)O[C@@H]%10[C@@H](COP(=O)(O)O[C@@H]%11[C@@H](COP(=O)(O)O[C@@H]%12[C@@H](COP(=O)(O)O[C@@H]%13[C@@H](COP(=O)(O)O[C@H]%14C[C@H](n%15cc(C)c(=O)[nH]c%15=O)O[C@@H]%14COP(=O)(O)O[C@@H]%14[C@@H](COP(=O)(O)O[C@@H]%15[C@@H](COP(=O)(O)O[C@H]%16C[C@H](n%17cc(C)c(=O)[nH]c%17=O)O[C@@H]%16COP(=O)(O)O[C@@H]%16[C@@H](COP(=O)(O)O[C@@H]%17[C@@H](COP(=O)(O)O[C@@H]%18[C@@H](COP(=O)(O)O[C@H]%19C[C@H](n%20cc(C)c(=O)[nH]c%20=O)O[C@@H]%19COP(=O)(O)O[C@@H]%19[C@@H](COP(=O)(O)O[C@@H]%20[C@@H](COc%21ccc%22c(c%21)Oc%21cc(O)ccc%21C%22%21OC(=O)c%22ccccc%22%21)O[C@@H](n%21cnc%22c(N)ncnc%22%21)[C@@H]%20O)O[C@@H](n%20cnc%21c(=O)[nH]c(N)nc%21%20)[C@@H]%19O)O[C@@H](n%19ccc(N)nc%19=O)[C@@H]%18O)O[C@@H](n%18ccc(N)nc%18=O)[C@@H]%17O)O[C@@H](n%17cnc%18c(=O)[nH]c(N)nc%18%17)[C@@H]%16O)O[C@@H](n%16cnc%17c(=O)[nH]c(N)nc%17%16)[C@@H]%15O)O[C@@H](n%15cnc%16c(=O)[nH]c(N)nc%16%15)[C@@H]%14O)O[C@@H](n%14cnc%15c(N)ncnc%15%14)[C@@H]%13O)O[C@@H](n%13cnc%14c(=O)[nH]c(N)nc%14%13)[C@@H]%12O)O[C@@H](n%12cnc%13c(=O)[nH]c(N)nc%13%12)[C@@H]%11O)O[C@@H](n%11cnc%12c(=O)[nH]c(N)nc%12%11)[C@@H]%10O)O[C@@H](n%10ccc(N)nc%10=O)[C@@H]9O)O[C@@H](n9cnc%10c(N)ncnc%109)[C@@H]8O)O[C@@H](n8cnc9c(=O)[nH]c(N)nc98)[C@@H]7O)O[C@@H](n7cnc8c(=O)[nH]c(N)nc87)[C@@H]6O)O[C@@H](n6cnc7c(=O)[nH]c(N)nc76)[C@@H]5O)O[C@@H](n5cnc6c(=O)[nH]c(N)nc65)[C@@H]4O)O[C@@H](n4cnc5c(=O)[nH]c(N)nc54)[C@@H]3O)O[C@@H](n3cnc4c(=O)[nH]c(N)nc43)[C@@H]2O)O1");
}
