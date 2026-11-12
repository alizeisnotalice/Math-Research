u[p_]:=InverseCDF[NormalDistribution[],1-p];
gaps=N[Table[u[p]-Sqrt[2] u[Sqrt[p]],{p,{1/4,1/2,3/4}}],40];
Print[gaps];If[!TrueQ[Min[gaps]>0],Exit[1]];
