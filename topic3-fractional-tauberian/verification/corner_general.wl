ClearAll[x,k,a]; f[a_]:=a^2/4+a/2+x(1-x)+2k^2-(k+Max[0,k-x+a/2])(k+Max[0,k-1+x+a/2]);
vals={f[-2x],f[2x],f[Clip[2(x-k),{-2x,2x}]],f[Clip[2(1-x-k),{-2x,2x}]]};
Do[Print[TimeConstrained[FullSimplify[Integrate[PiecewiseExpand[Max@@vals,0<x<1/2&&condition],{x,0,1/2},Assumptions->condition]*8,condition],40,"timeout"]],{condition,{0<k<1/2,1/2<k<1,k>1}}];
