m=1+ProductLog[1/E]; f[v_]=(Exp[v]+1)/v;
tests={FullSimplify[(D[f[v],v] v^2)==Exp[v](v-1)-1],
 FullSimplify[Integrate[Exp[-b(z-p)],{z,p-v/b,p}]+2Integrate[Exp[-b(z-p)],{z,p,1}]==(Exp[v]+1-2Exp[-b(1-p)])/b,b>0&&0<p<1&&v>0],
 Abs[N[f[m]-Exp[m],50]]<10^-40};
Print[tests];If[!And@@tests,Exit[1]];
