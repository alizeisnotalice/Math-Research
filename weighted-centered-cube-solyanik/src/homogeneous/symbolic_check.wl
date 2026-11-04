Print["Wolfram version: ",$Version];
Clear[y,z,n,tau,q];
i1=Integrate[(1+y^2)^(-1/2),{y,0,1}];
i2=Integrate[(1+y)^(-1),{y,0,1}];
i3=Integrate[(1+y+z)^(-2),{y,0,1},{z,0,1}];
checks=FullSimplify[{i1==ArcSinh[1],i2==Log[2],i3==Log[4/3]}];
Print["Face integrals: ",checks];
If[checks=!={True,True,True},Exit[1]];
factors={Sqrt[Pi]/(2 i1),1/(Sqrt[2] i2),6^(-2/3)/i3};
Print["Closed factors: ",factors];Print["Numerical factors: ",N[factors,35]];
window=FullSimplify[Equivalent[tau/(n q)<(n-2 tau)/(n (n-tau)),
 tau^2-(n+2 q)tau+n q>0],Assumptions->n>0&&q>0&&0<tau<n];
Print["General shift-order window: ",window];If[window=!=True,Exit[1]];
