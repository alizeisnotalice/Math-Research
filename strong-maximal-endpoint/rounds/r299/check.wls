Clear[x,s,q,h,j];
checks={FullSimplify[D[h[s Log[q]],{q,2}]==s^2/q^2 (h''[s Log[q]]-h'[s Log[q]]/s),Assumptions->s>0&&q>0],FullSimplify[D[x^j Exp[-2x],{x,2}]/(x^j Exp[-2x])==((j/x)-2)^2-j/x^2,Assumptions->x>0]};
Print[checks];If[!And@@checks,Exit[1]];
