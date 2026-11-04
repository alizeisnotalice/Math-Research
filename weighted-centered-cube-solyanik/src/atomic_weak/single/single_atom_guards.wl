Print["wolfram_version=",$Version];
Clear[t,r];
j=((r-1)^t+(r+1)^t)/(2 t r);n=(r-1)^t+(r+1)^t;
guards={
 TrueQ[FullSimplify[Limit[t*((t)^t+(2+t)^t)/(2*t*(1+t)),t->0,Direction->"FromAbove"]==1]],
 TrueQ[FullSimplify[D[j,r]==(r D[n,r]-n)/(2 t r^2)]],
 TrueQ[FullSimplify[r D[n,{r,2}]<0,Assumptions->0<t<1&&r>1]],
 TrueQ[FullSimplify[(((r-1)^(1/8)+(r+1)^(1/8))/(2*(1/8)*r)/.r->17/16)>27/4]]
 };
Print["guard_vector=",guards];If[!And@@guards,Exit[1]];Print["all_guards=True"];
