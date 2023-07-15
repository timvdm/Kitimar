#include "Validate.hpp"

void Rarey_smarts_part_3(OpenBabel::OBMol &mol)
{
    // SMARTS 101 - 150
    validate_contains<"C[OX2][CX3](=[OX1])[OX2]C">(mol);
    //validate_contains<"Cc1:c(O):c:c:[$(cCl),$([cH])]:c1">(mol); // FIXME: recursive
    //validate_contains<"Cl[CH2,CH]C(Cl)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    validate_contains<"C~C(~O)~O">(mol);
    validate_contains<"C~C~C~C">(mol);
    validate_contains<"C~C~O">(mol);
    //validate_contains<"F[CH2,CH]C(F)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    //validate_contains<"I[CH2,CH]C(I)[$(N(=O)~O),$(C=O),$(C#N),$(S=O),$(C(=O)N),$(a)]">(mol); // FIXME: recursive
    //validate_contains<"N!@[C$(C@*)]">(mol); // FIXME: recursive
    //validate_contains<"N!@[N$(N@*)]">(mol); // FIXME: recursive
    validate_contains<"N#CC(=O)">(mol);
    validate_contains<"N#CC[OH]">(mol);
    validate_contains<"N-!@N">(mol);
    validate_contains<"N-!@O">(mol);
    validate_contains<"N-!@P">(mol);
    validate_contains<"N-!@[Br]">(mol);
    validate_contains<"N-!@[Cl]">(mol);
    validate_contains<"N-!@[F]">(mol);
    validate_contains<"N-!@[I]">(mol);
    validate_contains<"N-@N">(mol);
    validate_contains<"N-@O">(mol);
    validate_contains<"N-@P">(mol);
    validate_contains<"N-@[Br]">(mol);
    validate_contains<"N-@[Cl]">(mol);
    validate_contains<"N-@[F]">(mol);
    validate_contains<"N-@[I]">(mol);
    validate_contains<"N-[C,c,N]=[C,c,N,n,O,S]">(mol);
    validate_contains<"N1CCC1=O">(mol);
    validate_contains<"N1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[O,N]">(mol);
    validate_contains<"N=!@N">(mol);
    validate_contains<"N=!@O">(mol);
    validate_contains<"N=!@P">(mol);
    validate_contains<"N=@N">(mol);
    validate_contains<"N=C=N">(mol);
    validate_contains<"N=C=[S,O]">(mol);
    validate_contains<"N=NC(S)N">(mol);
    validate_contains<"N=[O,C,N,S]">(mol);
    validate_contains<"NP(=O)(N)N">(mol);
    validate_contains<"N[CH2]C#N">(mol);
    validate_contains<"N[CX4H2][CX3](=[OX1])[O,N]">(mol);
    validate_contains<"N~C(~O)~N">(mol);
    validate_contains<"N~C~C(~O)~O">(mol);
    //validate_contains<"O!@[C$(C@*)]">(mol); // FIXME: recursive
    validate_contains<"O-!@O">(mol);
    validate_contains<"O-!@S">(mol);
    validate_contains<"O-!@[Br]">(mol);
    validate_contains<"O-!@[Cl]">(mol);
    validate_contains<"O-!@[F]">(mol);
    validate_contains<"O-!@[I]">(mol);
    validate_contains<"O-@O">(mol);
}
