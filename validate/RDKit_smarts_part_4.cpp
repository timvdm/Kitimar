#include "Validate.hpp"

void RDKit_smarts_part_4(OpenBabel::OBMol &mol)
{
    // SMARTS 151 - 200
    validate_contains<"[Si]">(mol);
    validate_contains<"[N;H2]-S">(mol);
    validate_contains<"[N;H]-[C;H]-[N;H]">(mol);
    validate_contains<"c-C(-O)(-O)-c">(mol);
    validate_contains<"c-C(-[O;H])[!O]">(mol);
    validate_contains<"c-S(=O)(=O)-C=C">(mol);
    validate_contains<"c-S(=O)(=O)-O">(mol);
    validate_contains<"n1(=O)c([F,Br,I,Cl])cccc1">(mol);
    validate_contains<"n1c([F,Br,I,Cl])cccc1">(mol);
    validate_contains<"n1c([F,Br,I,Cl])ncccc1">(mol);
    validate_contains<"C1-C-O-C-O-C1">(mol);
    validate_contains<"C=S">(mol);
    validate_contains<"c=S">(mol);
    validate_contains<"S=C-S">(mol);
    validate_contains<"C=C-N(=O)(=O)">(mol);
    validate_contains<"C=C-N(=O)(-O)">(mol);
    validate_contains<"C(=O)-C(=N)">(mol);
    validate_contains<"C=C-S">(mol);
    validate_contains<"C(-S)(-S)(-S)(-S)">(mol);
    validate_contains<"S(-O)(-O)(-O)(-O)">(mol);
    validate_contains<"S(=O)(=O)-C-N(=O)(=O)">(mol);
    validate_contains<"S(=O)(=O)-C-N(=O)(-O)">(mol);
    validate_contains<"N#C-S">(mol);
    validate_contains<"C=N-O">(mol);
    validate_contains<"C=N-S=N">(mol);
    validate_contains<"S-C(=O)-S">(mol);
    validate_contains<"S-C(=N)-S">(mol);
    validate_contains<"S-C(=N)(=N)">(mol);
    validate_contains<"N(=O)(=O)">(mol);
    validate_contains<"N(=O)(-O)">(mol);
    validate_contains<"N(-O)(-O)">(mol);
    validate_contains<"C1-S-C-S1">(mol);
    validate_contains<"N=C=S">(mol);
    validate_contains<"N=C=N">(mol);
    validate_contains<"C1-N=N1">(mol);
    validate_contains<"C([F,Br,Cl,I])-S">(mol);
    validate_contains<"C(-S)(-S)=C(-S)(-S)">(mol);
    validate_contains<"C(=N)-C=N-O">(mol);
    validate_contains<"C(-C#N)(-C#N)">(mol);
    validate_contains<"C-O-N=O">(mol);
    validate_contains<"C=C-N">(mol);
    validate_contains<"C=N-S(=O)(=O)">(mol);
    validate_contains<"N-N=O">(mol);
    validate_contains<"N-C(=O)-S">(mol);
    validate_contains<"S-C#N">(mol);
    validate_contains<"N-C(=O)-C(=O)-N">(mol);
    validate_contains<"c=N">(mol);
    validate_contains<"C1-O-C-O-C1">(mol);
    validate_contains<"C=N-N-S(=O)(=O)">(mol);
    validate_contains<"[C;!R](=O)-[C;!R](=O)">(mol);
}
