Print["Wolfram version: ",$Version];
phi[x_,p_]:=x^p-1-p(x-1);
Print["balanced norm expansion: ",FullSimplify[Series[(((1+d)^p+(1-d)^p)/2)^(2/p),{d,0,2}],Assumptions->p>1]];
Print["balanced S expansion: ",FullSimplify[Series[(((1+d)^p+(1-d)^p)/2)^(1/p),{d,0,2}],Assumptions->p>1]];
Do[
  varianceTopic22=q(1-q); holeDefectTopic22=q^(k+1)+(1-q)(1+q)^k-1;
  lim=Limit[holeDefectTopic22/varianceTopic22,q->0,Direction->"FromAbove"];
  Print["hole p=",k," limit D/var=",lim];If[!TrueQ[lim==k-1],Exit[1]];
  ep=M^(-k-1); c=(1-ep*M)/(1-ep); Ds=ep*M^k+(1-ep)c^k-1;
  distance=ep*(M-1)^k+(1-ep)*(1-c)^k;
  ls=Limit[distance/Ds,M->Infinity];Print["spike p=",k," limit=",ls];If[!TrueQ[ls==1],Exit[1]],{k,{2,3,4,8}}];
Print["p3 scalar factor:",Factor[phi[x,3]-2(x-1)^2]];
Print["p4 scalar factor:",Factor[phi[x,4]-3(x-1)^2]];
