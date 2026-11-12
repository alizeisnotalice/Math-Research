mu = 1 + ProductLog[1/E];
Do[t = N[InverseGammaRegularized[n, 1/8],24];
 If[Abs[GammaRegularized[n,t]-1/8]>10^-16,Exit[1]];
 Print[{n, N[(mu n-t+Log[1/4])/n,9]}],{n,{8,32,128,512,2048,8192}}];
