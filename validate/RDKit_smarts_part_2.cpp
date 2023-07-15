#include "Validate.hpp"

void RDKit_smarts_part_2(OpenBabel::OBMol &mol)
{
    // SMARTS 51 - 100
    validate_contains<"O=C1CCCC(N1)=O">(mol);
    validate_contains<"O1CCCCC1C2CCCO2">(mol);
    validate_contains<"[OH]c1cc([OH])cc2OC(C([OH])Cc21)c3cc([OH])c([OH])cc3">(mol);
    validate_contains<"C12OCCC(O1)CC2">(mol);
    validate_contains<"c-N(=O)~O">(mol);
    validate_contains<"[O,N,S][CH2]N1C(=O)CCC1(=O)">(mol);
    validate_contains<"O=C-N!@N">(mol);
    validate_contains<"C!@N=*">(mol);
    validate_contains<"[S,P](=O)(=O)OC">(mol);
    //validate_contains<"[$(N!@N),$([N;R0]=[N;R0])]">(mol); // FIXME: recursive
    //validate_contains<"[$(N!@S),$(N~C(~S)-N)]">(mol); // FIXME: recursive
    validate_contains<"[CH2]=[CH]-[N,O,S]">(mol);
    validate_contains<"S-C#N">(mol);
    validate_contains<"S=C-[#6,N,O]">(mol);
    validate_contains<"[Cl,Br,I]-N">(mol);
    //validate_contains<"O=CC([$(C(F)(F)F),$(C#N),$(Cl)])=C">(mol); // FIXME: recursive
    validate_contains<"C#C-[F,Br,I,Cl]">(mol);
    validate_contains<"C(-[O;H1])(-C#N)">(mol);
    validate_contains<"C(-O)-C-N(=O)=O">(mol);
    validate_contains<"C(=N)-S">(mol);
    validate_contains<"C(=O)-C(-O)-N-C=O">(mol);
    validate_contains<"C(=O)-C(=C)-C(=O)">(mol);
    validate_contains<"C(=O)-C(=N)">(mol);
    validate_contains<"C(=O)-C([F,Br,I,Cl])-C(=O)">(mol);
    validate_contains<"C(=O)-C([F,Br,I,Cl])=C([F,Br,I,Cl])-C(=O)">(mol);
    validate_contains<"C(=O)-C=C-C(=O)">(mol);
    validate_contains<"C(=O)-N(-O)-C(=O)">(mol);
    validate_contains<"C(=O)-N(-S)-C(=O)">(mol);
    validate_contains<"[C;!R](=O)-[N;!R]-[C;!R](=O)">(mol);
    validate_contains<"C(=O)-N-N(=O)">(mol);
    validate_contains<"C(=O)-N-N-C(=O)">(mol);
    validate_contains<"C(=O)-O-N-C(=O)">(mol);
    validate_contains<"C(=O)-S">(mol);
    validate_contains<"C(=O)-S-C(=S)">(mol);
    validate_contains<"C(=S)-S">(mol);
    validate_contains<"[C;H1](=O)">(mol);
    validate_contains<"C([F,Br,I,Cl])=N">(mol);
    validate_contains<"[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]-[C;H2]">(mol);
    validate_contains<"C-N=O">(mol);
    validate_contains<"C-S(=O)(=O)-O">(mol);
    validate_contains<"C1-C-C(=O)-O1">(mol);
    validate_contains<"C=C(-S)-S(=O)">(mol);
    validate_contains<"[C;!R]=[C;!R]-[C;!R](-O)">(mol);
    validate_contains<"C=C-C(=N)">(mol);
    validate_contains<"[C;!R]=[C;!R]-[C;!R](=O)">(mol);
    validate_contains<"C=C-C(=O)-C=C">(mol);
    validate_contains<"C=C-C(=S)">(mol);
    validate_contains<"C=C-C(=S)-S">(mol);
    validate_contains<"C=C-C([F,Br,I,Cl])=C([F,Br,I,Cl])">(mol);
    validate_contains<"C=C-C=N">(mol);
}
