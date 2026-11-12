m=1+ProductLog[1/E];u[p_]:=InverseCDF[NormalDistribution[],1-p];
d=CDF[NormalDistribution[],1/4]-1/2;
Print[N[{d,Sqrt[m]u[d]},18]];If[!TrueQ[0<N[d]<1/2],Exit[1]];
