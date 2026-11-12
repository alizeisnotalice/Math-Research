shared=Table[2Integrate[1-u,{u,a,a+1/2}],{a,{0,1/2}}];
prices=2-shared;Print[{shared,prices}];If[!TrueQ[prices=={5/4,7/4}],Exit[1]];
