r={1/3,1/3};tests={2(1/2)^2(2/3-1/3)^2==1/18,Total[r]^2-Total[r^2]==2/9,Total[r]-1/2==1/6};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
