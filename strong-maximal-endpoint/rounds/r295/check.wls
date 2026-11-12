Clear[s,b,d];
checks={FullSimplify[Integrate[Exp[-b/(12s)],{b,0,Infinity},Assumptions->s>0]==12s,Assumptions->s>0],Log[2]-3/2 < -1/4, FullSimplify[Sqrt[d]+d/Sqrt[d]==2 Sqrt[d],Assumptions->d>0]};
Print[checks];If[!And@@checks,Exit[1]];
