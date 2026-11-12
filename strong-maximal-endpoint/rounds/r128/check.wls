p=1/2; b=1/(3-y); bound=(1+b+1-p)/p;
tests={FullSimplify[1/3<b<p&&bound==3+2b,0<y<1],FullSimplify[2Integrate[bound-3,{y,0,1}]==4Log[3/2]]};
Print[{tests,N[2Log[3/2],20]}];If[!And@@tests,Exit[1]];
