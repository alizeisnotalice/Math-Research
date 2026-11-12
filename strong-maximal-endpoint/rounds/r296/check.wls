Clear[u,d,s,v,j];
checks={FullSimplify[Integrate[u (1+u^2/4)^(-d/4),{u,0,Infinity},Assumptions->d>4]==8/(d-4),Assumptions->d>4],
FullSimplify[Integrate[2u Exp[-u/s],{u,0,Infinity},Assumptions->s>0]==2s^2,Assumptions->s>0],
FullSimplify[D[j Log[v]-v,v,v]==-j/v^2,Assumptions->v>0]};
Print[checks];If[!And@@checks,Exit[1]];
