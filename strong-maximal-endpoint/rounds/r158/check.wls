cost[l_,h_]:=l-h+Binomial[h,2];
tests=Table[Min[Table[cost[l,h],{h,0,l}]]==l-1,{l,1,64}];Print[{And@@tests,1-1/16}];If[!TrueQ[And@@tests],Exit[1]];
