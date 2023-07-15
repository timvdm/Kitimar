#include "Validate.hpp"

void Rarey_smarts_part_24(OpenBabel::OBMol &mol)
{
    // SMARTS 1151 - 1200
    validate_contains<"c:1:c:c:c:c:c:1-[#7](-[#1])-[#6](=[#16])-[#7](-[#1])-[#6](-[#1])(-[#1])-c:2:n:c:c:c:c:2">(mol);
    //validate_contains<"c:2(:c:1-[#16]-c:3:c(-[#7](-c:1:c(:c(:c:2-[#1])-[#1])-[#1])-[$([#1]),$([#6](-[#1])(-[#1])-[#1]),$([#6](-[#1])(-[#1])-[#6]-[#1])]):c(:c(~[$([#1]),$([#6]:[#6])]):c(:c:3-[#1])-[$([#1]),$([#7](-[#1])-[#1]),$([#8]-[#6;X4])])~[$([#1]),$([#7](-[#1])-[#6;X4]),$([#6]:[#6])])-[#1]">(mol); // FIXME: recursive
    validate_contains<"c:2(:c:1-[#6](-[#6](-[#6](-c:1:c(:c(:c:2-[#1])-[#1])-[#1])(-[#1])-[#1])=[#8])=[#6](-[#6](-[#1])-[#1])-[#6](-[#1])-[#1])-[#1]">(mol);
    validate_contains<"c:2(:c:1:c(:c(:c(:c(:c:1:c(:c(:c:2-[#1])-[#1])-[#1])-[#1])-[#7](-[#1])-[#7](-[#1])-[#6]=[#8])-[#1])-[#1])-[#1]">(mol);
    validate_contains<"c:2(:c:1:c(:c(:c(:c(:c:1:c(:c(:c:2-[#8]-[#1])-[#6]=[#8])-[#1])-[#1])-[#1])-[#1])-[#1])-[#7](-[#1])-[#1]">(mol);
    validate_contains<"c:2(:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#6]=[#6]-[#6]-3=[#7])-[#7]">(mol);
    validate_contains<"c:2(:c:1:c:c:c:c:c:1:c-3:c(:c:2)-[#6](-c:4:c:c:c:c:c-3:4)=[#8])-[#8]-[#1]">(mol);
    validate_contains<"c:2(:c:1:c:c:c:c:c:1:n:n:c:2)-[#6](-[#6]:[#6])-[#6]#[#7]">(mol);
    validate_contains<"c:2(:n:c:1:c(:c(:c:c(:c:1-[#1])-[F,Cl,Br,I])-[#1]):n:2-[#1])-[#16]-[#6](-[#1])(-[#1])-[#6](=[#8])-[#7](-[#1])-[#6]:[#6]">(mol);
    validate_contains<"c:2-3:c(:c:c:1:c:c:c:c:c:1:c:2)-[#7](-[#6](-[#1])-[#1])-[#6](=[#8])-[#6](=[#7]-3)-[#6]:[#6]-[#7](-[#1])-[#6](-[#1])-[#1]">(mol);
    validate_contains<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7](-[#6;X4]-[#7]-3-[#1])-[#1]">(mol);
    validate_contains<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7](-[#7]=[#6]-3)-[#1]">(mol);
    validate_contains<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7]-[#6]=[#7]-3">(mol);
    validate_contains<"c:2:c:1:c:c:c:c-3:c:1:c(:c:c:2)-[#7]-[#7]=[#7]-3">(mol);
    validate_contains<"c:2:c:c:1:n:c(:c(:n:c:1:c:c:2)-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]:[#6])-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]:[#6]">(mol);
    validate_contains<"c:2:c:c:1:n:c(:c(:n:c:1:c:c:2)-c:3:c:c:c:c:c:3)-c:4:c:c:c:c:c:4-[#8]-[#1]">(mol);
    validate_contains<"c:2:c:c:1:n:c:3:c(:n:c:1:c:c:2):c:c:c:4:c:3:c:c:c:c:4">(mol);
    validate_contains<"c:2:c:c:1:n:n:c(:n:c:1:c:c:2)-[#6](-[#1])(-[#1])-[#6]=[#8]">(mol);
    validate_contains<"c:2:c:c:c:1:c(:c:c:c:1):c:c:2">(mol);
    validate_contains<"cC[N+]">(mol);
    validate_contains<"cN=[N+]=[N-]">(mol);
    validate_contains<"n(:c)(:c):a">(mol);
    validate_contains<"n-!@C">(mol);
    validate_contains<"n-!@N">(mol);
    validate_contains<"n-!@O">(mol);
    validate_contains<"n-!@S">(mol);
    validate_contains<"n-@C">(mol);
    validate_contains<"n-@N">(mol);
    validate_contains<"n-@O">(mol);
    validate_contains<"n-@P">(mol);
    validate_contains<"n-@S">(mol);
    validate_contains<"n-[Br]">(mol);
    validate_contains<"n-[Cl]">(mol);
    validate_contains<"n-[F]">(mol);
    validate_contains<"n-[I]">(mol);
    validate_contains<"n1(-[#6;X4])c(c(-[#1])c(c1-[#6]:[#6])-[#1])-[#6](-[#1])-[#1]">(mol);
    validate_contains<"n1(-[#6](-[#1])-[#1])c(c(-[#6](=[#8])-[#6])c(c1-[#6]:[#6])-[#6])-[#6](-[#1])-[#1]">(mol);
    validate_contains<"n1(-[#6])c(c(-[#1])c(c1-[#6](-[#1])(-[#1])-[#7](-[#1])-[#6](=[#16])-[#7]-[#1])-[#1])-[#1]">(mol);
    validate_contains<"n1(-[#6])c(c(-[#1])c(c1-[#6](-[#1])=[#6](-[#6]#[#7])-c:2:n:c:c:s:2)-[#1])-[#1]">(mol);
    validate_contains<"n1(-[#6])c(c(-[#1])c(c1-[#6](-[#1])=[#6]-2-[#6](=[#8])-[!#6&!#1]-[#6]=,:[!#1]-2)-[#1])-[#1]">(mol);
    validate_contains<"n1(-[#6])c(c(-[#1])c(c1-[#6]=[#7]-[#7])-[#1])-[#1]">(mol);
    validate_contains<"n1-2cccc1-[#6]=[#7](-[#6])-[#6]-[#6]-2">(mol);
    validate_contains<"n1nnnc2cccc12">(mol);
    validate_contains<"n1nscc1-c2nc(no2)-[#6]:[#6]">(mol);
    validate_contains<"n2(-[#6](-[#1])-[#1])c-1c(-[#6]:[#6]-[#6]-1=[#8])cc2-[#6](-[#1])-[#1]">(mol);
    validate_contains<"n2(-[#6](-[#1])-c:1:c(:c(:c:c(:c:1-[#1])-[#1])-[#1])-[#1])c(c(-[#1])c(c2-[#6]-[#1])-[#1])-[#6]-[#1]">(mol);
    validate_contains<"n2(-[#6]:1:[!#1]:[!#6&!#1]:[!#1]:[#6]:1-[#1])c(c(-[#1])c(c2-[#6;X4])-[#1])-[#6;X4]">(mol);
    validate_contains<"n2(-[#6]:1:[!#1]:[#6]:[#6]:[#6]:[#6]:1)c(cc(c2-[#6;X4])-[#1])-[#6;X4]">(mol);
    validate_contains<"n2(-[#6]:1:[!#1]:[#6]:[#6]:[#6]:[#6]:1)c(cc(c2-[#6]:[#6])-[#1])-[#6;X4]">(mol);
    validate_contains<"n2(-[#6]:1:[#6](-[#6]#[#7]):[#6]:[#6]:[!#6&!#1]:1)c(c(-[#1])c(c2)-[#1])-[#1]">(mol);
}
