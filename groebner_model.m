q := 11;
n := 17;
m := 251;
k := 16;
d:=Maximum((m-n),0);
//n1 := (k-1) * n;

F := GF(q);


// variables
if d eq 0 then
    R<[vars]> := PolynomialRing(F, 1*n + (k-1), "grevlex");
    lambda := [1] cat vars[1..(k-1)];
    kernel := Matrix(R, 1, m, vars[k..1*m + (k-1)]);
else 
    R<[vars]> := PolynomialRing(F, d*n + (k-1), "grevlex");
    lambda := [1] cat vars[1..(k-1)];
    kernel := HorizontalJoin(IdentityMatrix(R,d), Matrix(R, d, n, vars[k..d*n + (k-1)]));
end if;

// matrices
load "matrices.m";
Aone:=Matrix(F,m,n,ElementToSequence(A[1]));
AprimeOne, EchTransform := EchelonForm(Aone);
AprimeOne:=Matrix(R,m,n,ElementToSequence(AprimeOne));
EchTransform:=Matrix(R,m,m,ElementToSequence(EchTransform));
Aprime := [EchTransform * A[i] - lambda[i] * AprimeOne : i in [2..k]]; 
Ajoined := HorizontalJoin(Aprime);


system := ElementToSequence(kernel[1]*Ajoined);  
// system cat:= // fill in 


SetVerbose("Groebner", 1);
GB := GroebnerBasis(system);
GB;
exit;


