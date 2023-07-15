#include "Validate.hpp"

void BindingDB_substructure_part_54(OpenBabel::OBMol &mol)
{
    // SMARTS 2651 - 2700
    //validate_contains<"[H]N(N([H])C1=NC(=NC(=C1[H])C([H])([H])[H])N([H])S(=O)(=O)C2=C([H])C([H])=C([H])C([H])=C2[H])C([H])=C3/C([H])=C([H])C(=O)C(OC([H])([H])[H])=C3[H]">(mol); // FIXME: stereo
    validate_contains<"[H]N(N([H])S(=O)(=O)C1=C([H])C([H])=C(N([H])C(=O)C([H])([H])[H])C([H])=C1[H])C(=O)C([H])([H])OC2=C(Cl)C([H])=C(Cl)C([H])=C2[H]">(mol);
    validate_contains<"[H]N([H])C">(mol);
    //validate_contains<"[H]N([H])C(=N/C1=NC2=C([H])C([H])=C([H])C([H])=C2S1)N([H])C(=O)C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H]">(mol); // FIXME: stereo
    validate_contains<"[H]N([H])C(=O)C1=C(C=C([H])C([H])=C1)N([H])C2=NC(=NC([H])=C2F)N([H])C3=C([H])C(=C(C([H])=C3[H])N4C([H])([H])C([H])([H])N(C([H])([H])[H])C([H])([H])C4([H])[H])C([H])([H])[H]">(mol);
    validate_contains<"[H]N([H])C(=O)C1=C(C=C([H])C([H])=C1)N([H])C2=NC(=NC([H])=C2F)N([H])C3=C([H])C(=C(N)C([H])=C3[H])C([H])([H])[H]">(mol);
    validate_contains<"[H]N([H])C(=O)C1=CC=CC=C1">(mol);
    validate_contains<"[H]N([H])C(C)=O">(mol);
    validate_contains<"[H]N([H])C1=CC=CC=C1">(mol);
    validate_contains<"[H]N([H])C1=NC(=O)C(=N1)C2=CC([H])([H])N([H])C(=O)C3=C2C([H])=C(Br)N3[H]">(mol);
    validate_contains<"[H]N([H])C1=NC2=C(N=C([H])N2[H])C(O)=N1">(mol);
    validate_contains<"[H]N([H])C1=NC=CC=*1C">(mol);
    validate_contains<"[H]N([H])C1=NC=CC=C1">(mol);
    validate_contains<"[H]N([H])C1=NC=CC=C1C">(mol);
    validate_contains<"[H]N([H])C1=NCCN1C">(mol);
    validate_contains<"[H]N([H])C1=NCCO1">(mol);
    validate_contains<"[H]N([H])CC1CCCC1">(mol);
    validate_contains<"[H]N([H])S(=O)(=O)C1=C([H])C([H])=C(N([H])N([H])C2=C(O)N([H])C3=C2C([H])=C(Br)C([H])=C3[H])C([H])=C1[H]">(mol);
    validate_contains<"[H]N([H])S(=O)(=O)C1=C([H])C([H])=C(N([H])N([H])C2C(=O)N([H])C3=C2C([H])=C(Br)C([H])=C3[H])C([H])=C1[H]">(mol);
    //validate_contains<"[H]N([H])S(=O)(=O)C1=C([H])C([H])=C(N([H])N=C2/C(=O)N([H])C3=C2C([H])=C(Br)C([H])=C3[H])C([H])=C1[H]">(mol); // FIXME: stereo
    validate_contains<"[H]N([H])S(=O)(=O)C1=CC=C(NC2=NC([H])=C([H])C=N2)C=C1">(mol);
    validate_contains<"[H]N([H])S(=O)(=O)C1=NC([H])=C(CC2=CC(C)=CC3=C2C(=O)N(C3Cl)C4=CC=C(F)C(O)=C4)C([H])=C1[H]">(mol);
    //validate_contains<"[H]N([H])[C@@]1([H])C([H])([H])C([H])([H])[C@]([H])(N([H])C2=NC(N([H])C([H])([H])C3=C([H])C([H])=C([H])C([H])=C3[H])=C4N=C([H])N(C4=N2)C5([H])C([H])([H])C([H])([H])C([H])([H])C5([H])[H])C([H])([H])C1([H])[H]">(mol); // FIXME: stereo
    validate_contains<"[H]N1C(=O)CCCCC1([H])[H]">(mol);
    validate_contains<"[H]N1C(Br)=C([H])C2=C1C(=O)N([H])C([H])([H])C=C2">(mol);
    validate_contains<"[H]N1C([H])=C2C(=C1O)C3=C(NC4=C3C([H])=C([H])C([H])=C4[H])C5=C2C6=C(N5)C([H])=C([H])C([H])=C6[H]">(mol);
    validate_contains<"[H]N1C([H])=CC([H])=C([H])C1=O">(mol);
    validate_contains<"[H]N1C=C(F)C(=O)N([H])C1=O">(mol);
    validate_contains<"[H]N1CC2=C(CC1C(O)=O)C=CC=C2">(mol);
    validate_contains<"[H]N1CC2=CC=CC=C2C3=CC=CN=C13">(mol);
    validate_contains<"[H]N1CCN(CC1)C2=CC=C(NC3=NC4=C(C=N3)C([H])=C([H])C(=O)N4)*=C2">(mol);
    validate_contains<"[H]N1N=C(C(CC2=CC=CC=C2)=C1C(=O)NC)C3=CC=CC=C3">(mol);
    validate_contains<"[H]N1N=C(C)C(CC)=C1C(=O)NC">(mol);
    validate_contains<"[H]N=C(N([H])C)N([H])C">(mol);
    validate_contains<"[H]NC(=O)N([H])C(C)=O">(mol);
    validate_contains<"[H]NC1=CC=C(O[H])C=C1">(mol);
    validate_contains<"[H]NS(=O)(=O)C1=CC2=CC=CC=C2C=C1">(mol);
    validate_contains<"[H]OC(=O)C(=O)C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]">(mol);
    validate_contains<"[H]OC(=O)C1=CC=CC=C1OC(C)=O">(mol);
    validate_contains<"[H]OC(=O)CC(NC(C)=O)C(N)=O">(mol);
    //validate_contains<"[H]OC(=O)[C@@]1([H])O[C@@]([H])(OC2=C3C(=O)C([H])=C(OC3=C(C(O[H])=C2[H])[C@@]4([H])C([H])([H])C([H])([H])N(C([H])([H])[H])C([H])([H])[C@@]4([H])O[H])C5=C([H])C([H])=C([H])C([H])=C5Cl)[C@]([H])(O[H])[C@@]([H])(O[H])[C@]1([H])O[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]OC([H])([H])C(=O)[C@@]1([H])C([H])([H])C([H])([H])[C@@]2([H])[C@]3([H])C([H])([H])C([H])([H])C4=C([H])C(=O)C([H])([H])C([H])([H])[C@]4(C([H])([H])[H])[C@@]3([H])[C@@]5([H])O[C@]([H])(O[H])[C@@]12C5([H])[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]OC([H])([H])[C@]([H])(n1c([H])nc(c1[H])C(=O)N([H])[H])C([H])([H])C([H])([H])n2c([H])c([H])c3c([H])c([H])c(N([H])C(=O)C([H])([H])C([H])([H])c4c([H])c([H])c([H])c([H])c4[H])c([H])c23">(mol); // FIXME: stereo
    validate_contains<"[H]OC1=C(O[H])C2=C3C(=C1[H])C(=O)OC1=C3C(=C([H])C(O[H])=C1O[H])C(=O)O2">(mol);
    //validate_contains<"[H]OC1=C([H])C([H])=C([H])C(C(=O)N([H])[C@@]([H])(C([H])([H])C2=C([H])C([H])=C([H])C([H])=C2[H])[C@]([H])(O[H])C(=O)N3C([H])([H])SC(C([H])([H])[H])(C([H])([H])[H])[C@@]3([H])C(=O)N([H])C([H])([H])C4=C([H])C([H])=C([H])C([H])=C4C([H])([H])[H])=C1C([H])([H])[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]OC1=C([H])C2=C(C([H])=C1[H])[C@@]3([H])C([H])([H])C([H])([H])[C@@]4(C(=O)C([H])([H])C([H])([H])[C@@]4([H])[C@]3([H])C([H])([H])C2([H])[H])C([H])([H])[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]OC1=C([H])C2=C(C([H])=C1[H])[C@@]3([H])C([H])([H])C([H])([H])[C@]4(C([H])([H])[H])[C@@]([H])(O[H])C([H])([H])C([H])([H])[C@@]4([H])[C@]3([H])C([H])([H])C2([H])[H]">(mol); // FIXME: stereo
    validate_contains<"[H]OC1=CC=C(Cl)C(O[H])=C1">(mol);
    validate_contains<"[H]OC1=CC=CC2=C1C(CCN(C)C)=CN2">(mol);
    validate_contains<"[H]OC1C([H])C([H])C2C(C1[H]3CC(C4CCC(CC4)C5CCCCC5)C6OC7C(CCC(C8CCCCC8)C7OC6C3C9CCCCC9)CCCC(CC)CCCCCC)C([H])([H])C([H])([H])C([H])C2([H])C([H])([H])C([H])[H]C([H])([H])CC(O)C([H])([H])C([H])([H])C[H]">(mol);
}
