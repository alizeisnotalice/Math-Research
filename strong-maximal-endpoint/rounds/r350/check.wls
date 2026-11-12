Clear[h,l,n];checks={FullSimplify[Integrate[Exp[-l h],{h,0,Infinity},Assumptions->l>0]==1/l,Assumptions->l>0],Limit[n^2 ((1+1/n^3)^(n-1)-1),n->Infinity]==1};Print[checks];If[!And@@checks,Exit[1]];
