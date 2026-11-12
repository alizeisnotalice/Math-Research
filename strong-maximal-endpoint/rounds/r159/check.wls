r={9/16,1/16};k=Total[r]^2-Total[r^2];
checks={Total[r]==5/8,k==9/128,k>=Total[r](Total[r]-Max[r])};Print[checks];If[!TrueQ[And@@checks],Exit[1]];
