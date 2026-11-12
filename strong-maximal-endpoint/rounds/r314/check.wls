Clear[v,t,j,x];
checks=Table[FullSimplify[D[v Exp[-(t+Log[v])] Sum[(t+Log[v])^j/j!,{j,0,k-1}],v]==Exp[-(t+Log[v])] Sum[(t+Log[v])^j/j!,{j,0,k-2}],Assumptions->v>0],{k,2,15}];
Print[{Length[checks],And@@checks}];If[!And@@checks,Exit[1]];
