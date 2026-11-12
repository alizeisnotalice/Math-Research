ratios=Join[Table[2^i,{i,-5,16}],Table[3*2^i/2,{i,-5,16}]];
tests=Table[Count[Range[0,20],j_/;2^j<x]==Ceiling[Max[0,Log[2,x]]],{x,ratios}];
Print[{Length[tests],And@@tests}];If[!TrueQ[And@@tests],Exit[1]];
