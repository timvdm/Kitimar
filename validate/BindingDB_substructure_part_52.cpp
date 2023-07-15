#include "Validate.hpp"

void BindingDB_substructure_part_52(OpenBabel::OBMol &mol)
{
    // SMARTS 2551 - 2600
    validate_contains<"S=C1NC=CN1">(mol);
    validate_contains<"S=C1NC=NN1">(mol);
    validate_contains<"S=CCC=S">(mol);
    validate_contains<"S=NC1=CC=CC=C1">(mol);
    validate_contains<"SC#N">(mol);
    validate_contains<"SC1=CC2=C(C=CC=C2)C=C1">(mol);
    validate_contains<"SC1=CC=CC=C1Br">(mol);
    validate_contains<"SC1=CNC2=C1C=CC=C2">(mol);
    validate_contains<"SC1=NC2=CC=CC=C2N1">(mol);
    validate_contains<"SC1=NC=CC=C1Cl">(mol);
    validate_contains<"SC1=NN=C(N1C2=CC=CC=C2)C3=CC=NC=C3">(mol);
    validate_contains<"SC1=NN=CN1C2=CC=CC=C2">(mol);
    validate_contains<"SC=N">(mol);
    validate_contains<"SCCCC(=O)NC1CN(CC2=CC=CC=C2)CCN1">(mol);
    validate_contains<"[C++](C1=CC=CC=C1)C2=CC=CC=C2">(mol);
    validate_contains<"[C+3]C([C+3])(O)O">(mol);
    //validate_contains<"[CH3:26][C@H:12]([C:10](=[CH2:11])[CH2:9][C:3]1=[CH:5][CH2:4][C:1](=[CH:2]1)[CH:6]2[CH2:8][CH2:7]2)[C:13]3=[CH:16][CH:15]=[C:14]([CH:18]=[CH:17]3)[CH:22]4[CH2:20][CH2:21][C@@H:24]([CH3:25])[C:19]4=[CH2:23]">(mol); // FIXME: stereo
    validate_contains<"[CH3:29][O:28][C:15]1=[C:10]([C:8](=[O:9])[C:3]2=[C:4]([NH2:7])[N:5]=[C:1]([NH:16][CH:17]3[CH2:18][CH2:19][N:22]([CH2:20][CH2:21]3)[S:24]([CH3:23])(=[O:25])=[O:26])[N:6]=[CH:2]2)[C:11]([F:27])=[C:12]([F:30])[CH:13]=[CH:14]1">(mol);
    //validate_contains<"[CH3:37][O:7][C:29]1=[CH:27][CH:21]=[C:16]([CH:22]=[CH:28]1)[C@@:9]23[O:1][C:15]4=[C:12]([C:17]([O:4][CH3:33])=[CH:23][C:24]([O:6][CH3:36])=[CH:20]4)[C@:10]2([OH:2])[C@H:13]([OH:3])[C@@H:14]([C@H:11]3[C:18]5=[CH:25][CH:30]=[CH:32][CH:31]=[CH:26]5)[C:19](=[O:5])[N:8]([CH3:34])[CH3:35]">(mol); // FIXME: stereo
    //validate_contains<"[ClH:1].C1=CC=[C+:1]C=C1">(mol); // FIXME: components
    //validate_contains<"[FH2+].CC.CN.CC1=CC=C(C=C1)[S]([O-])[O-]">(mol); // FIXME: components
    //validate_contains<"[H+:2].[H+:2].[OH2:2]">(mol); // FIXME: components
    //validate_contains<"[H+].C.N">(mol); // FIXME: components
    //validate_contains<"[H+].C1=CC2=C(C=C1)N=CC=C2">(mol); // FIXME: components
    //validate_contains<"[H+].CC12CCC3C(CCC4=C3C=CC(O)=C4)C1CCO2">(mol); // FIXME: components
    //validate_contains<"[H+].N.O">(mol); // FIXME: components
    //validate_contains<"[H+].O.O.OO.OO.OC1C[H]C[H]1">(mol); // FIXME: components
    //validate_contains<"[H+].[H+].O.CON(CC(C)C)SC1=CC=NC=C1">(mol); // FIXME: components
    //validate_contains<"[H+].[H+].[H+].N">(mol); // FIXME: components
    //validate_contains<"[H:31][N:6]([H:32])[C@@:8]1([H:24])[C:10]([H:26])([H:27])[C:12](=[C:11]([H:28])[C@@:9]([H:25])([O:1][C:13]([H:29])([C:15]([H:35])([H:36])[C:19]([H:40])([H:41])[H:42])[C:14]([H:33])([H:34])[C:18]([H:37])([H:38])[H:39])[C@:7]1([H:23])[N:5]([H:30])[C:17](=[O:4])[C:20]([H:43])([H:44])[H:45])[C:16](=[O:3])[O:2][C:21]([H:46])([H:47])[C:22]([H:48])([H:49])[H:50]">(mol); // FIXME: stereo
    //validate_contains<"[H:788][N:63]([H:789])[C:60](=[O:788])[C:57](=[O:788])[C:54](=[O:787])[N:50]([H:787])[C:22]([H:786])([H:786])[C:21]([H:46])([H:47])[O:2][C:16](=[O:3])[C:12]1=[C:11]([H:28])[C@@:9]([H:25])([O:1][C:13]([H:29])([C:15]([H:35])([H:36])[C:19]([H:40])([H:41])[H:42])[C:14]([H:33])([H:34])[C:18]([H:37])([H:38])[H:39])[C@:7]([H:23])([N:5]([H:30])[C:17](=[O:4])[C:20]([H:43])([H:44])[H:45])[C@@:8]([H:24])([N:6]([H:31])[H:32])[C:10]1([H:26])[H:27]">(mol); // FIXME: stereo
    validate_contains<"[H]C">(mol);
    validate_contains<"[H]C(=C(C)C1=CC=C(O)C=C1)C2=CC(O)=CC(O)=C2">(mol);
    validate_contains<"[H]C(=C)C1=C([H])NC2=C1C([H])=CC=C2[H]">(mol);
    validate_contains<"[H]C(=O)C(O)COP(O)(O)=O">(mol);
    validate_contains<"[H]C(=O)C1=C(O)C(O)=CC(Br)=C1Br">(mol);
    validate_contains<"[H]C(=O)C1=CC2=C(OCCO2)C=C1">(mol);
    validate_contains<"[H]C(=O)N(C)O">(mol);
    validate_contains<"[H]C(=O)OC1=CC=CN=N1">(mol);
    //validate_contains<"[H]C(=O)[C@H](CCCNC(N)=N)NC(=O)[C@H](CCCN)NC(=O)[C@H](CCCCN)NC(=O)CC1=CC=CC=C1">(mol); // FIXME: stereo
    validate_contains<"[H]C(C)(Br)C([H])(C)Br">(mol);
    validate_contains<"[H]C(C)=NO">(mol);
    validate_contains<"[H]C(CCCC(=O)N([H])C1=C(N)C=CC=C1)CC(=O)C2=CC3=CC=CC=C3C=C2">(mol);
    validate_contains<"[H]C(Cl)=C=C([H])Cl">(mol);
    validate_contains<"[H]C([H])(N)C([H])([H])N1C=CN=C1">(mol);
    validate_contains<"[H]C([H])(N)C([H])([H])S(N)(=O)=O">(mol);
    validate_contains<"[H]C([H])(N)CC(CC(CC1=CC=C(O)C=C1)CC([H])([H])CC([H])([H])CC(CC(CC(C)C)C(O)=O)C(C)CC)CC(O)=O">(mol);
    validate_contains<"[H]C([H])(S)C([H])(C)C(N)=O">(mol);
    validate_contains<"[H]C([H])([H])C(N)C(O)=O">(mol);
    //validate_contains<"[H]C([H])([H])[H].[H]OC(=O)[C@@]([H])(N([H])C(=O)C([H])([H])[H])C([H])([H])C1=C([H])N([H])C2=C1C([H])=C([H])C([H])=C2[H]">(mol); // FIXME: components
}
