v = 3 + 2 Integrate[4/(1+x)-1,{x,0,1}];
If[!TrueQ[Simplify[v==1+8 Log[2]]],Exit[1]];
Print[{v,N[3/v,20]}]
