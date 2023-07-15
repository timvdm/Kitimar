#include "Validate.hpp"

void BindingDB_substructure_part_55(OpenBabel::OBMol &mol)
{
    // SMARTS 2701 - 2750
    validate_contains<"[H]OO[H]">(mol);
    //validate_contains<"[H]OP(=O)(O[H])O[C@]([H])(C([H])([H])N(C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H])S(=O)(=O)C1=C([H])C([H])=C(N([H])[H])C([H])=C1[H])[C@@]([H])(N([H])C(=O)O[C@]2([H])C([H])([H])OC([H])([H])C2([H])[H])C([H])([H])C3=C([H])C([H])=C([H])C([H])=C3[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]O[C@]([H])(C([H])([H])N(C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H])S(=O)(=O)C1=C([H])C([H])=C(C([H])=C1[H])N([H])([H])=N)[C@@]([H])(N([H])C(=O)O[C@]2([H])C([H])([H])OC([H])([H])C2([H])[H])C([H])([H])C3=C([H])C([H])=C([H])C([H])=C3[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]O[C@]([H])(C([H])([H])N(C([H])([H])C([H])(C([H])([H])[H])C([H])([H])[H])S(=O)(=O)C1=C([H])C([H])=C(N([H])[H])C([H])=C1[H])[C@@]([H])(N([H])C(=O)O[C@]2([H])C([H])([H])OC([H])([H])C2([H])[H])C([H])([H])C3=C([H])C([H])=C([H])C([H])=C3[H]">(mol); // FIXME: stereo
    //validate_contains<"[H]O[C@]1([H])c2c([H])c3OC([H])([H])Oc3c([H])c2[C@]([H])(c4c([H])c(OC([H])([H])[H])c(OC([H])([H])[H])c(OC([H])([H])[H])c4[H])[C@]5([H])C(=O)OC([H])([H])[C@]15[H]">(mol); // FIXME: stereo
    validate_contains<"[H]SC([H])([H])C([H])(N([H])[H])C([H])([H])C([H])([H])N([H])[H]">(mol);
    //validate_contains<"[H][C@@](CCCC[NH+](C)C)(NC(=O)[C@@]([H])(NC(=O)[C@]([H])(CC1=CC=C(C=C1)C(N)=[NH2+])NC(C)=O)C2CCCCC2)C(=O)N[C@@]([H])(CC(C)C)C(=O)N3CCC[C@@]3([H])C(N)=O">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@](O)(COP(O)(O)=O)C=O">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CC3=CNC4=CC=CC(=C34)[C@@]1([H])C[C@H](CN2CC=C)C(N)=O">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CCC(C)=C[C@@]1([H])C3=C(O)C=C(CCCCC)C=C3OC2(C)C">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CCC3CCCC(=O)[C@]3([H])C1C[C@H](OC2=O)[C@@H]4CCOC4">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CCCC[C@]1(C)C1CC[C@]3(C)CCCC3C1CC2">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CCCC[C@]1([H])CN(C[C@@H](O)[C@H](CSC3=CC=CC=C3)NC(=O)C4=C(C)C(O)=CC=C4)[C@@H](C2)C(=O)NC(C)(C)C">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]3([H])[C@@]2([H])CCC4=CC(=O)CC[C@]34C">(mol); // FIXME: stereo
    //validate_contains<"[H][C@@]12NC(=O)N[C@]1([H])[C@H](CCCCC(O)=O)SS2">(mol); // FIXME: stereo
    //validate_contains<"[H][C@](CC1=CC=CC=C1)(NC(=O)OCC)C(=O)N2CCC[C@@]2([H])C(=O)NNCCCC[NH3+]">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]1(CCN(CC2=CC=CC=C2)C1)NS(=O)(=O)NCCC3=C(NC(=N3)C4=CC=CC=C4)C5=CC=C(OC)C=C5">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]1(COP(O)(=O)OP([O-])(=O)OP(O)([O-])=O)O[C@@]([H])(N2C=CC(N)=NC2=O)C([H])(O)[C@@]1([H])O">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]1(COP(O)(=O)OP([O-])(=O)OP(O)([O-])=O)O[C@@]([H])(N2C=NC3=C2N=CN=C3N)C([H])(O)[C@@]1([H])O">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]1(N[C@@H](C(=O)NCCNC(=O)[C@@H]2N[C@]([H])(SC2(C)C)[C@H](NC(=O)CC3=CC=C4C=CC=CC4=C3)C(=O)NCC)C(C)(C)S1)[C@H](NC(=O)CC5=CC6=C(C=CC=C6)C=C5)C(=O)NCC">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]12CC[C@]([H])(C1)[C@@H](C2)C(=O)NC3=C(C(=O)OC)C(CC)=C(C)S3">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]12CS[C@@H](CCCCC(O)=O)[C@@]1([H])NC(=O)N2">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]12OCC[C@@]1([H])[C@H](CO2)OC(=O)N[C@@H](Cc1ccc(OCc2ccccn2)cc1)[C@H](O)CN(CC(C)C)S(=O)(=O)c1ccc2OCOc2c1">(mol); // FIXME: stereo
    //validate_contains<"[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@H](N)C3=CC=CC=C3)C(O)=O">(mol); // FIXME: stereo
    validate_contains<"[H][N+]1(C)C=C(CCN)C2=C1NC3=C(CCC3C(N)O)C2=O">(mol);
    validate_contains<"[H][N+]1=C(C=C(C=C1NCCC2=CC=CC=C2)C3=CC=CC=C3)C4=CC=C(O)C=C4">(mol);
    validate_contains<"[H][N+]1=C(N)C2=C(NC(=N2)[S](CCCCC)C3=NC4=C(S3)C(Cl)=CC=C4)N=C1">(mol);
    //validate_contains<"[H]\C">(mol); // FIXME: stereo
    validate_contains<"[NH-]C1=NC=CC=N1">(mol);
    validate_contains<"[NH2:1]CN">(mol);
    //validate_contains<"[NH2:7][C:2]1=[N:6][C:3]([NH2:11])=[C:4]([N:8]=[O:9])[C:1]([O:10][CH2:12][C@H:13]2[CH2:18][CH2:17][CH:16]=[CH:15][CH2:14]2)=[N:5]1">(mol); // FIXME: stereo
    validate_contains<"[NH3+]C(CCCCB(O)O)C([O-])=O">(mol);
    //validate_contains<"[NH3+][C@H](CC1=CC=CC=C1)C(=O)N2CCC[C@H]2C(=O)NCC3=CC=CC(F)=C3">(mol); // FIXME: stereo
    //validate_contains<"[NH][C@@H](CCCCN)C=O">(mol); // FIXME: stereo
    validate_contains<"[Na]Cl">(mol);
    validate_contains<"[O-]C1=CC=CC=C1">(mol);
    validate_contains<"[O-]O[S](N1CCNCC1)C2=CC=CC3=C2C=CC=C3Cl">(mol);
    validate_contains<"[O-]S(=O)(=O)CC[NH+]1CCOCC1">(mol);
    //validate_contains<"[O-]S(=O)(=O)c1ccc2NC(=O)C(c2c1)=C3/Nc4ccccc4C3=O">(mol); // FIXME: stereo
    validate_contains<"[O-][N+](=O)C1=CC=CN1">(mol);
    validate_contains<"[O-][N+](=O)C1=CC=CO1">(mol);
    validate_contains<"[O-][N+](=O)C1=CC=CS1">(mol);
    //validate_contains<"[OH-].[OH-].CC12CCCC1C3CCC4=CC=CC=C4C3CC2">(mol); // FIXME: components
    //validate_contains<"[W].[H]C([H])([H])C([H])([H])[H].[H]C([H])([H])C([H])([H])[H].[H]C([H])([H])C([H])([H])[H].[H]C1([H])[W]C([H])([H])C([H])([H])C1([H])[H].[H]C2([H])C([H])([H])C([H])([H])C([H])([H])C2([H])[H].[H]C3([H])C([H])([H])C([H])([H])C([H])([H])C3([H])[H].[H]C([H])([H])C4([H])C([H])([H])C([H])([H])C([H])([H])C4([H])[H]">(mol); // FIXME: components
    validate_contains<"[Zn]">(mol);
    validate_contains<"[nH]1cnc4c1c2c(nc(cc2Nc3ccc(cc3)N(C)C(=O)C)C)cc4">(mol);
    validate_contains<"c1cc(N)ccc1">(mol);
    validate_contains<"c1cc(N)ccc1N">(mol);
    validate_contains<"c1cc(N)ccc1OS">(mol);
    validate_contains<"c1cc2cc[nH]c2[nH]1">(mol);
}
