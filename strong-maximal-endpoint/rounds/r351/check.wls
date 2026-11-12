ClearAll[tail, quantile];
mu = 1 + ProductLog[1/E];
q = (mu - 1)/mu;
tail[n_, d_, h_?NumericQ] := Sum[
  PDF[BinomialDistribution[n, N[q,32]], j] GammaRegularized[n+d, Max[0,h-j N[mu,32]]], {j,0,n}];
quantile[n_, d_, p_] := h /. FindRoot[tail[n,d,h]==p,{h,N[n mu,24]},WorkingPrecision->24];
Do[a=quantile[n,0,p]; b=quantile[n,1,p];
 If[!(b>a && Abs[tail[n,1,b]-p]<10^-15),Exit[1]];
 Print[{n,p,N[b-a,10]}],{n,{8,32,128}},{p,{1/4,1/2,3/4}}];
