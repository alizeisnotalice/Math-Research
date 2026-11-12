f[g_]=Integrate[2(4/(1+z)-1),{z,g/2,g},Assumptions->0<g<1];
h[g_]=Integrate[2(4/(1+z)-1),{z,g/2,1},Assumptions->1<g<2];
ok={FullSimplify[D[f[g],g]>0,0<g<1],FullSimplify[D[h[g],g]<0,1<g<2],FullSimplify[f[1]==h[1]==8Log[4/3]-1]};
Print[{ok,N[f[1]/2,20]}];If[!And@@ok,Exit[1]];
