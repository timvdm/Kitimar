#include "Validate.hpp"

void Rarey_smarts_part_1(OpenBabel::OBMol &mol)
{
    // SMARTS 1 - 50
    validate_contains<"*!:[a]:*:*:[a]!:*">(mol);
    validate_contains<"*!@*">(mol);
    validate_contains<"*!@[#8]!@*">(mol);
    validate_contains<"*(!@*)(!@*)">(mol);
    validate_contains<"*-!:[a]:*:*:[a]-!:*">(mol);
    validate_contains<"*-!:[a]:*:[a]-!:*">(mol);
    validate_contains<"*-!:aa-!:*">(mol);
    validate_contains<"*-!:aaa-!:*">(mol);
    validate_contains<"*-!:aaaa-!:*">(mol);
    validate_contains<"*-!@*">(mol);
    validate_contains<"*-*-[R]:[a]">(mol);
    validate_contains<"*-*-[R]~*:[a]">(mol);
    validate_contains<"*-*-[R]~*~*:[a]">(mol);
    //validate_contains<"*/,\[R]=,:;@[R]/,\*">(mol); // FIXME: stereo
    //validate_contains<"*/,\[R]=;@[R]/,\*">(mol); // FIXME: stereo
    validate_contains<"*1*******1">(mol);
    validate_contains<"*1******1">(mol);
    validate_contains<"*1*****1">(mol);
    validate_contains<"*1****1">(mol);
    validate_contains<"*1***1">(mol);
    validate_contains<"*1**1">(mol);
    validate_contains<"A-[Cl,Br,I]">(mol);
    //validate_contains<"A[#6;!$(C(F)(F)F)]">(mol); // FIXME: recursive
    //validate_contains<"Br[CH2,CH]C(Br)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    //validate_contains<"C!@[C$(C@*)]">(mol); // FIXME: recursive
    //validate_contains<"C!@[N$(C@*)]">(mol); // FIXME: recursive
    //validate_contains<"C!@[N$(N@*)]">(mol); // FIXME: recursive
    //validate_contains<"C!@[N$(O@*)]">(mol); // FIXME: recursive
    //validate_contains<"C!@[N$(S@*)]">(mol); // FIXME: recursive
    validate_contains<"C#!@C">(mol);
    validate_contains<"C#!@N">(mol);
    validate_contains<"C#@C">(mol);
    //validate_contains<"C#C[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    validate_contains<"C(=*)(=*)">(mol);
    validate_contains<"C(=O)(~c)~c">(mol);
    validate_contains<"C(=O)-N-C=O">(mol);
    validate_contains<"C(=O)C(=O)">(mol);
    validate_contains<"C(=O)C[N+,n+]">(mol);
    validate_contains<"C(=O)N(C(=O))OC(=O)">(mol);
    validate_contains<"C(=O)OC(=O)">(mol);
    validate_contains<"C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)">(mol);
    validate_contains<"C(=O)Oc1ccc(N(=O)=O)cc1">(mol);
    validate_contains<"C(=O)Onnn">(mol);
    validate_contains<"C(=O)[Cl,Br,I]">(mol);
    validate_contains<"C(=O)[OH1]">(mol);
    validate_contains<"C(=O)[SH1]">(mol);
    validate_contains<"C(=[O,S])(N)Oc">(mol);
    validate_contains<"C(C)(C)=CC=[O,SX2]">(mol);
    //validate_contains<"C([F,Cl,Br,I])[CH2,CH][$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    validate_contains<"C-!@C">(mol);
}
