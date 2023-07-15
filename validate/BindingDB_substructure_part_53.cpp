#include "Validate.hpp"

void BindingDB_substructure_part_53(OpenBabel::OBMol &mol)
{
    // SMARTS 2601 - 2650
    validate_contains<"[H]C([H])C([H])([H])S(N)(=O)=O">(mol);
    validate_contains<"[H]C([H])C2(CCC1=CC=C(N)C=C1)CC(O)=C(SC3=C(C=C(CO)C([H])=C3)C([H])([H])[H])C(=O)O2">(mol);
    validate_contains<"[H]C1(O)NC=NC2=C1N=CN2C3OC(O)C(O)C3C">(mol);
    validate_contains<"[H]C1(OC(CO)C(O)C(O)C1O)C2=CC=C(Cl)C(CC3=CC=C(OCC)C=C3)=C2">(mol);
    //validate_contains<"[H]C1([H])[C@@]([H])(O)[C@@]([H])(COP(O)(=O)OP(O)(=O)OP(O)(=O)O[Na])O[C@@]1([H])N1C=C(C)C(=O)NC1=O">(mol); // FIXME: stereo
    validate_contains<"[H]C12CCOC1([H])OCC2OC(=O)NC(Cc3ccccc3)C(O)CN(CC(C)C)S(=O)(=O)c4ccccc4">(mol);
    validate_contains<"[H]C1=C(NC)N=C(NC)N=C1">(mol);
    validate_contains<"[H]C1=C(NC2=CC=CC=C2)C3=CC=CC=C3N=C1C">(mol);
    validate_contains<"[H]C1=C([H])C(=NC(NC)=N1)C2=C([H])C([H])=C(C)C([H])=C2[H]">(mol);
    validate_contains<"[H]C1=C([H])C(N)=NC(N)=N1">(mol);
    validate_contains<"[H]C1=C([H])C(O)=C(C([H])=C1[H])C(O)=O">(mol);
    validate_contains<"[H]C1=C([H])C2=C(NC1=O)N=C(N)N=C2">(mol);
    validate_contains<"[H]C1=C([H])C2=C(NC1=O)N=C(NC3=CC=CC=C3)N=C2">(mol);
    validate_contains<"[H]C1=C([H])C2=CN=C(N)N=C2NC1=O">(mol);
    validate_contains<"[H]C1=C([H])N=C(S1)C(=O)C([H])([H])[H]">(mol);
    validate_contains<"[H]C1=CC(N)=NC(N)=N1">(mol);
    validate_contains<"[H]C1=CC([H])=C([H])C(=O)N1">(mol);
    validate_contains<"[H]C1=CC2=NC([H])=C([H])C([H])=C2C([H])=C1[H]">(mol);
    validate_contains<"[H]C1=CN=C(NC)N=C1NC">(mol);
    validate_contains<"[H]C1=CN=CN=C1[H]">(mol);
    validate_contains<"[H]C1=NC(C)=NC(N)=C1[H]">(mol);
    validate_contains<"[H]C1=NC(N)=NC(C)=C1[H]">(mol);
    validate_contains<"[H]C1=NC(N)=NC(N)=C1[H]">(mol);
    validate_contains<"[H]C1=NC(NC)=NC(=C1[H])N(C)C">(mol);
    validate_contains<"[H]C1=NC(NC)=NC(NC)=C1[H]">(mol);
    validate_contains<"[H]C1=NC2=C(C=C(NC(=O)CC)C=C2)C(NC3=CC=CC=C3)=C1[H]">(mol);
    validate_contains<"[H]C1=NC2=C(C=CC=C2)C(NC)=C1[H]">(mol);
    validate_contains<"[H]C1=NC2=CC(O)=C(O)C=C2C(NC3=CC=CC=C3)=N1">(mol);
    validate_contains<"[H]C1=NC2=CC=CC=C2C(NC3=CC=CC=C3)=C1">(mol);
    //validate_contains<"[H]C1C([H])[C@@]([H])(OC([H])([H])[H])[C@@]2(O[C@@]1([H])N3C4=C(C([H])=C([H])C([H])=C4[H])C5=C3C6=C(C7=C([H])N([H])C(O)=C57)C8=C(C([H])=C([H])C([H])=C8[H])N26)C([H])([H])[H]">(mol); // FIXME: stereo
    validate_contains<"[H]C1CC(O)C2(C)CCC3C4=CC=C(O)C=C4CCC3([H])C12[H]">(mol);
    validate_contains<"[H]C1N([H])C(=O)C2C1C3=C(NC4=C3C([H])=C([H])C([H])=C4[H])C5=C2C6=C(N5)C([H])=C([H])C([H])=C6[H]">(mol);
    validate_contains<"[H]C1NC([H])C(C([H])C1[H])C(O)=O">(mol);
    validate_contains<"[H]CN(C)C(=O)C1=CC=CC=C1">(mol);
    validate_contains<"[H]COC1=C(O)C=CC(=C1)C([H])=O">(mol);
    validate_contains<"[H]N(C(=O)C(C(C)C)N(C)C(C)=O)C(C)(C)COC(=O)OC">(mol);
    validate_contains<"[H]N(C(=O)C1=CC=C2N=C(N)SC2=C1)C3=CC=CC=C3">(mol);
    validate_contains<"[H]N(C(C)=O)C1=CC(C)=C(C=C1)N(C)C">(mol);
    validate_contains<"[H]N(C(C)=O)C1=CC=CC=C1">(mol);
    validate_contains<"[H]N(C(C)=O)C1=CC=CC=C1Cl">(mol);
    validate_contains<"[H]N(C(C)=O)C1=CC=CC=C1N">(mol);
    validate_contains<"[H]N(C([H])([H])C([H])([H])[H])C([H])(C([H])([H])[H])C([H])([H])C1=C([H])C(=C([H])C([H])=C1[H])C(F)(F)F">(mol);
    validate_contains<"[H]N(C)C1CCN(C1)C1=CC=CC=C1">(mol);
    validate_contains<"[H]N(C1=CC=C(C=C1OC2=CC=CC=C2)[N+]([O-])=O)S(C)(=O)=O">(mol);
    validate_contains<"[H]N(C1=CC=CC=C1)C1=CC=CC2=C1C(=CC=C2)S(O)(=O)=O">(mol);
    validate_contains<"[H]N(C1=CC=CC=C1)C2=NC=NC3=CC=CC=C23">(mol);
    validate_contains<"[H]N(C1=NC([H])=NC2=C1N=C([H])N2C([H])([H])[H])C3=C([H])C(=C([H])C([H])=C3[H])C([H])([H])[H]">(mol);
    validate_contains<"[H]N(C1=NC=CC=N1)C2=C(C)C=CC=C2">(mol);
    validate_contains<"[H]N(CC(=O)NC(CC1=CC=CC=C1)C(=O)N2CCCCC2)S(=O)(=O)C3=CC4=CC=CC=C4C=C3">(mol);
    validate_contains<"[H]N(CC)CC">(mol);
}
