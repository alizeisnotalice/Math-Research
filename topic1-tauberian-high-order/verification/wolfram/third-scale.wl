Print["VERSION ",$Version];
f=(1+2u)^(-n);Print["BINOMIAL ",InputForm[FullSimplify[Table[D[f,{u,k}]/.u->0,{k,0,2}],Assumptions->n>=2 && Element[n,Integers]]]];
Print["q^2 order after /u^(n+2): ",InputForm[2^(2n) u^(n-2)]];
Do[Module[{q=(2u/(1+2u))^n,expr,coeff,expect},expr=Normal[Series[q/(1-q),{u,0,n+2}]];coeff=Coefficient[expr,u,n+2];expect=2^(n+1)n(n+1)+If[n==2,16,0];Print[InputForm[{"n"->n,"series"->expr,"third"->coeff,"expected"->expect,"check"->(coeff===expect)}]]],{n,2,10}];
