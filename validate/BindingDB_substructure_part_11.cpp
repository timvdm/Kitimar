#include "Validate.hpp"

void BindingDB_substructure_part_11(OpenBabel::OBMol &mol)
{
    // SMARTS 501 - 550
    validate_contains<"OC1CCCCC1C1=CC=C(C#N)C(=C1)C(F)(F)F">(mol);
    validate_contains<"OC1COC(OOP(O)(=O)OP(O)(O)=O)C1O">(mol);
    validate_contains<"OC1NC2=C(C=C1)C(O)=CC=C2">(mol);
    validate_contains<"OC1NC2=C(C=C1)C(OC1=CC=CC=C1)=CC=C2">(mol);
    validate_contains<"OCCN1C(=O)C2=CC=CC=C2N=C1C=CC1=CC(O)=C(O)C=C1">(mol);
    validate_contains<"OCCN=C">(mol);
    validate_contains<"ON1C(=O)C(=CC2=CN=C(NC3=CC=CC=C3)N=C12)C1=C(Cl)C=CC=C1Cl">(mol);
    validate_contains<"ON1C(=O)C(=CC2=CN=CN=C12)C1=CC=CC=C1">(mol);
    //validate_contains<"ON=C1/C(=O)NC2=CC=C(O)C=C12">(mol); // FIXME: stereo
    validate_contains<"ONC(=O)C([*])N([*])S(=O)(=O)C1=CC=C([*])C=C1">(mol);
    //validate_contains<"O[C@@H]1CC[C@]2(O)C[C@H]1OC2=O">(mol); // FIXME: stereo
    //validate_contains<"O[C@@H]1C[C@@](O)(C[C@@H](O)[C@H]1O)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"O[SH](O)=O">(mol);
    validate_contains<"Oc1cc(O)c(cc1Cl)-c1[nH]ncc1-c1ccc(OCCCC#N)cc1">(mol);
    validate_contains<"S1C2=C(C=CC=C2)N=C1C1=CC=CC=C1">(mol);
    validate_contains<"S1C=CC2=C1C=CC=N2">(mol);
    validate_contains<"S1C=CC2=C1C=CN=C2">(mol);
    validate_contains<"S1C=CC2=C1C=NC=C2">(mol);
    validate_contains<"S1C=CN2C=CC=CC12">(mol);
    validate_contains<"S1C=CN2C=CN=CC12">(mol);
    validate_contains<"S1C=NC(=C1)N1C=CC2=CC=CC=C12">(mol);
    validate_contains<"S1C=NC(=C1)N1C=NC2=CC=CC=C12">(mol);
    validate_contains<"S1C=NC(=C1)N1C=NC=N1">(mol);
    validate_contains<"S1C=NC=C1N1C=NC2=CC=CC=C12">(mol);
    validate_contains<"[*]N1C([*])=NC2=CC3=NC(=CN=C3C=C12)C1=CC=CC=C1">(mol);
    validate_contains<"[CH3:31][O:30][c:27]1[cH:26][cH:25][c:24]([cH:29][cH:28]1)-[c:23]1[o:32][c:33]2[n:34][cH:39][n:38][c:36]([NH2:37])[c:35]2[c:22]1-[c:19]1[cH:18][cH:17][c:16]([NH:15][C:13](=[O:14])[NH:12][c:11]2[cH:6][c:5]([cH:7][cH:8][c:9]2[F:10])[C:2]([F:1])([F:3])[F:4])[cH:21][cH:20]1">(mol);
    validate_contains<"[H:43][N:14]([c:15]1[n:16][c:17]([H:44])[c:18]([o:19]1)-[c:20]1[c:21]([H:45])[c:22]([H:46])[c:23]([H:47])[c:24]([c:25]1[H:48])-[c:26]1[n:31][c:30]([H:52])[c:29]([H:51])[c:28]([H:50])[c:27]1[H:49])[c:12]1[c:13]([H:42])[c:6]([c:7]([H:37])[c:8]([H:38])[c:9]1[O:10][C:11]([H:39])([H:40])[H:41])[S:2](=[O:1])(=[O:3])[C:4]([H:32])([H:33])[C:5]([H:34])([H:35])[H:36]">(mol);
    //validate_contains<"[H]C(=C(/[H])C(=O)OC([H])(C([H])([H])[H])C([H])([H])[H])N1N=C(N=C1[H])C1=C([H])C([H])=C(C([H])=C1[H])C1=C([H])C([H])=C([H])C([H])=C1C(F)(F)F">(mol); // FIXME: stereo
    validate_contains<"[H]C1(CCCCC1)OC">(mol);
    validate_contains<"[H]C1([H])CCC([H])([H])[N+]1([H])[H]">(mol);
    //validate_contains<"[H]C1=C([H])C([H])=C(OC([H])([H])[C@@]2([H])N(C([H])([H])[H])C([H])([H])C([H])([H])C2([H])[H])C(=N1)C([H])([H])[H]">(mol); // FIXME: stereo
    validate_contains<"[H]C1CCC([H])([H])[N+]1([H])[H]">(mol);
    validate_contains<"[H]N(C)C1=NN([H])C2=CC=CC=C12">(mol);
    validate_contains<"[H]N([H])C1=C([H])C([H])=C([H])C([H])=C1S(=O)(=O)N([H])[H]">(mol);
    validate_contains<"[H]N([H])C1CCNC1">(mol);
    validate_contains<"[H]N([H])CC1CC(C)C(C)C1">(mol);
    validate_contains<"[H]N([H])CC1CC(C)CC1C">(mol);
    validate_contains<"[H]N1CCN(C)CC1">(mol);
    //validate_contains<"[H]OC([H])([H])[C@]([H])(N([H])C(=O)[C@@]([H])([N+]([H])([H])[H])C([H])([H])C1=C([H])C([H])=C([H])C([H])=C1[H])C(=O)N([H])[C@]([H])(C([O-])=O)C([H])([H])C([H])([H])C([H])([H])N([H])C(N([H])[H])=[N+]([H])[H]">(mol); // FIXME: stereo
    validate_contains<"[H]OC1=C(C([H])=C([H])C([H])=C1[H])C(=O)OC([H])([H])[H]">(mol);
    //validate_contains<"[H][C@@]12CC[C@H](O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])C(CCCCCCCC(=O)NC)C[C@@]2([H])CC(=O)CC[C@]12C">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]1(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O[C@@]([H])(N2C=NC3=C2N=C(N)NC3=O)C([H])(O)[C@@]1([H])O">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]1(CO[P@@](O)(=O)O[P@@]([O-])(=O)O[P@](O)([O-])=O)O[C@@]([H])(N2C=NC3=C2N=CN=C3N)[C@@]([H])(O)[C@@]1([H])O">(mol); // FIXME: stereo
    validate_contains<"[H][N+]1([H])CCCC1">(mol);
    validate_contains<"[O-]C(=O)CN1C=CC2=CC=CC=C12">(mol);
    validate_contains<"[O-][N+](=O)C1=CC=C(C=C1)C1=NC(=C(N1)C1=CC=NC=C1)C1=CC=C(F)C=C1">(mol);
    validate_contains<"[O-][N+](=O)C1=CC=C(C=C1)C1=NC(=C(N1)C1=CC=NC=C1)C1=CC=CC=C1">(mol);
    validate_contains<"[O-][N+](=O)C1=CC=CC(CN2C(=O)C3=CC=CC=C3C3=CC=CC=C23)=C1">(mol);
    validate_contains<"[O-][S+]1(=O)NNC(=O)c2ccccc12">(mol);
    validate_contains<"*C(=O)NC1CC1">(mol);
}
