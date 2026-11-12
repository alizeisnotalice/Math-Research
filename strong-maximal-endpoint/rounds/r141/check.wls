cost[l_,r_]:=With[{a=Quotient[l,r],b=Mod[l,r]},(r-b)a^2+b(a+1)^2];
tests=Flatten[Table[cost[l,r]==Min[Total[#^2]& /@ Select[IntegerPartitions[l],Length[#]<=r&]],{l,1,16},{r,1,l}]];
Print[{Length[tests],And@@tests,cost[13,5]/13}];If[!TrueQ[And@@tests],Exit[1]];
