Clear[s,x,b,u];
moment=Integrate[(2/(1+x))^s,{x,1,Infinity},Assumptions->s>1];
checks={FullSimplify[moment==2/(s-1),Assumptions->s>1],Integrate[Exp[-b/6],{b,0,Infinity}]==6,
Sum[(2/3)^u,{u,0,Infinity}]==3};
Print[checks];If[!And@@checks,Exit[1]];
