Clear[x,s,l];
checks={FullSimplify[D[(s-1) Log[1+x]-(s-1)x,x]/.x->0]==0,
FullSimplify[D[x-Log[1+x],{x,2}]==1/(1+x)^2,Assumptions->x> -1],
FullSimplify[Integrate[Exp[-(s-1)x^2],{x,-Infinity,Infinity},Assumptions->s>1]==Sqrt[Pi/(s-1)],Assumptions->s>1]};
Print[checks];If[!And@@checks,Exit[1]];
